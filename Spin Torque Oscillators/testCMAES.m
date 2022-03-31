%function xmin=purecmaes   % (mu/mu_w, lambda)-CMA-ES
close all 
clear
rng(1,'twister');

config.res_type = 'multiMM';            % state type of reservoir(s) to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network with multiple functions), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more. Place reservoirs in cell ({}) for heterotic systems.
config.num_nodes = {[49]};                   % num of nodes in each sub-reservoir, e.g. if config.num_nodes = [10,5,15], there would be 3 sub-reservoirs with 10, 5 and 15 nodes each.

config = selectReservoirType(config);         % collect function pointers for the selected reservoir type

%% Evolutionary parameters
config.num_tests = 1;                         % num of tests/runs
config.test = 1;

config.error_to_check = 'train&val&test';

%% Task parameters
config.dataset = 'narma_10';          % Task to evolve for
config.figure_array = [figure figure];

% get any additional params. This might include:
% details on reservoir structure, extra task variables, etc.
config = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

% --------------------  Initialization --------------------------------
% User defined input parameters (need to be edited)

N = 3 +config.num_nodes{1}*2; %variables to manipulate [input scaling, W scaling , leak rate]
xmean = [2*rand-1; rand; rand; zeros(config.num_nodes{1}*2,1)];    % objective variables initial point

sigma = 0.3;          % coordinate wise standard deviation (step size)
stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
stopeval = 1e3*N^2;   % stop after stopeval number of function evaluations

% Strategy parameter setting: Selection
lambda = 4;%+floor(3*log(N));  % population size, offspring number
config.pop_size = lambda;
mu = lambda/2;               % number of parents/points for recombination
weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
mu = floor(mu);
weights = weights/sum(weights);     % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4+mueff/N) / (N+4 + 2*mueff/N);  % time constant for cumulation for C
cs = (mueff+2) / (N+mueff+5);  % t-const for cumulation for sigma control
c1 = 2 / ((N+1.3)^2+mueff);    % learning rate for rank-one update of C
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff));  % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
% usually close to 1

% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
B = eye(N,N);                       % B defines the coordinate system
D = ones(N,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of
%   ||N(0,I)|| == norm(randn(N,1))

% -------------------- Generation Loop --------------------------------
population = config.createFcn(config);
% make all the same
for i = 1:length(population)
    population(i) = population(1);
    population(i).core_indx = i;
end

minimum_weight = 0.1;
gen = 1;
counteval = 0;  % the next 40 lines contain the 20 lines of interesting code
best_error = inf;
while counteval < stopeval
    
    % Generate and evaluate lambda offspring
    parfor k=1:lambda
        
        tmp_arx = xmean + sigma * B * (D .* randn(N,1)); % m + sig * Normal(0,C)
        
        % evaluate function
        population(k).layer(1).input_scaling =  tmp_arx(1);
        population(k).layer(1).damping =  tmp_arx(2);
        population(k).layer(1).leak_rate =  tmp_arx(3);
        
        input_weights = reshape(tmp_arx(4:end),2,config.num_nodes{1});
        input_weights(input_weights<minimum_weight & input_weights~=0 & input_weights>-minimum_weight) = 0;
        population(k).layer(1).input_weights{1} = input_weights;

        population(k) = config.testFcn(population(k),config);
        arfitness(k) = getError(config.error_to_check,population(k));% objective function call
        
        fprintf('      Pop %d, Error: %.4f \n',k,arfitness(k));
        counteval = counteval+1;
        
        arx(:,k) = tmp_arx;
    end
    
    
    
    % Sort by fitness and compute weighted mean into xmean
    [arfitness, arindex] = sort(arfitness); % minimization
    xold = xmean;
    xmean = arx(:,arindex(1:mu))*weights;   % recombination, new mean value
    
    % Cumulation: Update evolution paths
    ps = (1-cs)*ps ...
        + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4 + 2/(N+1);
    pc = (1-cc)*pc ...
        + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
    
    % Adapt covariance matrix C
    artmp = (1/sigma) * (arx(:,arindex(1:mu))-repmat(xold,1,mu));
    C = (1-c1-cmu) * C ...                  % regard old matrix
        + c1 * (pc*pc' ...                 % plus rank one update
        + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
        + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
    
    % Adapt step size sigma
    sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
    
    % Decomposition of C into B*diag(D.^2)*B' (diagonalization)
    if counteval - eigeneval > lambda/(c1+cmu)/N/10  % to achieve O(N^2)
        eigeneval = counteval;
        C = triu(C) + triu(C,1)'; % enforce symmetry
        [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
        D = sqrt(diag(D));        % D is a vector of standard deviations now
        invsqrtC = B * diag(D.^-1) * B';
    end
    
    store_ardx(gen,:) = arx(:, arindex(1));
    store_err(gen,:) = population(arindex(1)).test_error;
    
    set(0,'currentFigure',config.figure_array(1))
    subplot(1,2,1)
    plot(store_ardx)
    subplot(1,2,2)
    plot(store_err)
    drawnow
    
    if store_err(gen,:) < best_error
        best_indv = population(arindex(1));
        plotReservoirDetails(best_indv,1,1,1,config)
        best_error = store_err(gen,:);
    end
    
    fprintf('Gen %d, Error: %.4f, Best Error: %.4f \n',gen,store_err(gen,:), min(min(store_err)));
            
    % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable
    if arfitness(1) <= stopfitness || max(D) > 1e7 * min(D)
        break;
    end
    
    gen = gen + 1;
end % while, end generation loop

xmin = arx(:, arindex(1)); % Return best point of last iteration.

function individual = updateParams(individual,config)

% end of first subres
end_pos = individual.nodes(1);

% intialise W and W_scale
individual.W_scaling(1:end_pos,1:end_pos) = individual.init_W_scaling(1);
individual.input_scaling(:,1:end_pos) = individual.init_input_scaling(:,1);

if config.mulit_leak_rate % assign leak rates for each node, if used
    individual.leak_rate = individual.init_leak_rate;
else
    individual.leak_rate(1:end_pos) = individual.init_leak_rate(1);
end

% update subres scaling
for i = 2:length(individual.nodes)
    % set subres positions
    start_pos = end_pos+1;
    end_pos = start_pos + individual.nodes(i)-1;
    
    if ~config.mulit_leak_rate
        individual.input_scaling(:,start_pos:end_pos) = individual.init_input_scaling(:,i);
    end
    individual.W_scaling(start_pos:end_pos,start_pos:end_pos) = individual.init_W_scaling(i);
    individual.leak_rate(start_pos:end_pos) = individual.init_leak_rate(i);
end

% indentify any W that are new connections and put W_scaling as identity
individual.W_scaling((individual.W.*individual.test_mask{1}) ~= 0) = 1;
% check zero connections do not have a scaling
individual.W_scaling(logical(individual.test_mask{1}.*((individual.W.*individual.test_mask{1}) == 0))) = 0;

end