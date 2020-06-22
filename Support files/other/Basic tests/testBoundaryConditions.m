clear
close all

%set random seed for experiments
rng(1,'twister');

%% Setup
config.parallel = 1;                        % use parallel toolbox

%start paralllel pool if empty
if isempty(gcp) && config.parallel %
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

% type of network to evolve
config.res_type = 'MM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = [100];                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

% Network details
config.metrics = {'KR','GR','linearMC'};       % behaviours that will be used; name metrics to use and order of metrics
config.voxel_size = 10;                  % when measuring quality, this will determine the voxel size. Depends on systems being compared. Rule of thumb: around 10 is good

% dummy variables for dataset; not used but still needed for functions to
% work
config.train_input_sequence= [];
config.train_output_sequence =[];
config.dataset = 'blank';

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = 10;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
comb = [0 0 0;
    0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 1 0;
    1 1 1];

config.test = 1;

% create population of reservoirs
population = config.createFcn(config);

ppm = ParforProgMon('Initial population: ', config.pop_size*length(comb));
for c = 1:length(comb)
    t_config = config;
    t_config.periodic_boundary = comb(c,:);
    
    for pop_indx = 1:config.pop_size
        warning('off','all')
        tic;
        t_config.periodic_boundary
        behaviours(c,pop_indx,:) = getMetrics(population(pop_indx),t_config);
        time_taken(c,pop_indx) = toc;
        ppm.increment();
    end
end

% thickness vs time
% figure
% plot(depth,time_taken)
% xlabel('Thickness (nm)')
% ylabel('Seconds')

% metrics vs thickness
figure
base_value = reshape(behaviours(1,:,:),10,3);

diff = behaviours-behaviours(1,:,:);

hold on
for b = 1:length(comb)
    plot(reshape(behaviours(b,:,:),10,3)')
end
hold off
legend('KR','GR','MC')
xlabel('Boundary Condition')
ylabel('Metric Value')