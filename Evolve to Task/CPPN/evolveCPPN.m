%% Evolve CPPN to construct substrate for a specific task
% This script can be used to evolve any reservoir directly to a task. It
% uses the steady-state Microbial Genetic Algorithm to evolve the best
% solution.

% Author: M. Dale
% Date: 08/11/18
clear
close all

% add all subfolders to the path --> make all functions in subdirectories available
% addpath(genpath(pwd));

warning('off','all')
rng(1,'twister');

%% Setup
config.parallel = 0;                        % use parallel toolbox

%start paralllel pool if empty
if isempty(gcp) && config.parallel
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

%% Evolutionary parameters
config.num_tests = 1;                        % num of runs
config.pop_size = 100;                       % large pop better
config.total_gens = 2000;                    % num of gens
config.mut_rate = 0.05;                       % mutation rate
config.deme_percent = 0.1;                  % speciation percentage
config.deme = round(config.pop_size*config.deme_percent);
config.rec_rate = 0.5;                       % recombination rate

%% substrate details
config_sub = config;
config_sub.pop_size = config.pop_size;
config_sub.num_nodes = [100];                 % num of nodes in subreservoirs, e.g. config.num_nodes = {10,5,15}, would be 3 subreservoirs with n-nodes each
config_sub.res_type ='RoRmin';               % currently only works with lattice as CPPN configured substrate
config_sub = selectReservoirType(config_sub);       % get correct functions for type of reservoir

%% CPPN details
config.res_type = 'RoRmin';                      % can use different hierarchical reservoirs. RoR_IA is default ESN.
config.num_nodes = [10];                  % num of nodes in subreservoirs, e.g. config.num_nodes = {10,5,15}, would be 3 subreservoirs with n-nodes each
config = selectReservoirType(config);       % get correct functions for type of reservoir
config.CPPN_inputs = 3;                     % coord of node A(x,y) and node B(x,y)
config.CPPN_outputs = length(config_sub.num_nodes)*2 +1; % Output 1: input layer, 2: hidden layer, 3: outputlayer
config.CPPN_list = {'internal'};

%config.preprocess = 1;                   % basic preprocessing, e.g. scaling and mean variance
config.dataset = 'CPPN';                 % Task to evolve for
% get any additional params
config = getAdditionalParameters(config);
config = selectDataset(config);

% add CPPN params
config.sparse_input_weights = 0;
config.sparsity = 1; 
config.add_input_states = 0;
config.multi_activ = 1;                      % use different activation funcs
config.activ_list = {@sin,@cos,@tan,@tanh,@biSigmoid,@gauss,@ramp,@step,@spike,@Inverse}; 
config.evolve_output_weights = 0;             % evolve rather than train
config.internal_weight_initialisation = 'norm';  % e.g.,  'norm', 'uniform', 'orth', etc.  must be same length as number of subreservoirs
config.evolve_output_weights = 1;
config.output_connectivity = 1;
%config.output_weight_scaler = 100;              % defines maximum/minimum weight value when evolving output weights

%% Task parameters
config_sub.dataset = 'narma_10';                 % Task to evolve for
% get dataset
config_sub = getAdditionalParameters(config_sub);
config_sub = selectDataset(config_sub);
% turn off certain parameters
config_sub.add_input_states = 1;
config_sub.CPPN_on = 1; % turn off weight mutation

%% general params
config.gen_print = 5;                       % gens to display achive and database
config.start_time = datestr(now, 'HH:MM:SS');
config.figure_array = [figure figure];
config_sub.figure_array = [figure figure];
config.save_gen = inf;                      % save at gen = saveGen
config.multi_offspring = 0;                  % multiple tournament selection and offspring in one cycle
config.num_sync_offspring = config.deme;      % length of cycle/synchronisation step
config.metrics = {'KR','GR','linearMC'};          % metrics to use
config.record_metrics = 0;                  % save metrics
config.error_to_check = 'train&val&test';

%% RUn MicroGA
for test = 1:config.num_tests
    
    clearvars -except config config_sub test storeError 
    
    fprintf('\n Test: %d  ',test);
    fprintf('Processing genotype......... %s \n',datestr(now, 'HH:MM:SS'))
    tic
    
    rng(test,'twister');
    
    last_best = 1;
    
    % create initial population
    CPPN = config.createFcn(config);
    
    % create initial substrate
    substrate = config_sub.createFcn(config_sub);
    
    %Assess population
    if config.parallel
        ppm = ParforProgMon('Initial population: ', config.pop_size);
        parfor pop_indx = 1:config.pop_size
            warning('off','all')
            % assign weights through CPPN
            [substrate(pop_indx),~,CPPN(pop_indx),~] = assessCPPNonSubstrate(substrate(pop_indx),config_sub,CPPN(pop_indx),config);
            ppm.increment();
            %fprintf('\n i = %d, error = %.4f, took: %.4f\n',popEval,substrate(popEval).valError,toc);
        end
    else
        for pop_indx = 1:config.pop_size
            % assign weights through CPPN
            tic
            [substrate(pop_indx),~,CPPN(pop_indx),~] = assessCPPNonSubstrate(substrate(pop_indx),config_sub,CPPN(pop_indx),config);
            fprintf('\n i = %d, error = %.4f, took: %.4f\n',pop_indx,getError(config.error_to_check,substrate(pop_indx)),toc);  
        end
    end
        
    % find an d print best individual
    [best(1),best_indv(1)] = min(getError(config.error_to_check,substrate));
    fprintf('\n Starting loop... Best error = %.4f\n',best);
    
    % store error that will be used as fitness in the GA
    store_error(test,1,:) = getError(config.error_to_check,substrate);
    
    %% start GA
    for gen = 2:config.total_gens
        
        % define seed
        rng(gen,'twister');
        
        % reshape stored error to compare
        cmp_error = reshape(store_error(test,gen-1,:),1,size(store_error,3));
        
        % Tournment selection - pick two individuals
        equal = 1;
        while(equal)
            indv1 = randi([1 config.pop_size]);
            indv2 = indv1+randi([1 config.deme]);
            if indv2 > config.pop_size
                indv2 = indv2- config.pop_size;
            end
            if indv1 ~= indv2
                equal = 0;
            end
        end
        
        % Assess fitness of both and assign winner/loser - highest score
        % wins
        if cmp_error(indv1) < cmp_error(indv2)
            winner=indv1; loser = indv2;
        else
            winner=indv2; loser = indv1;
        end
        
        % Manipulate CPPN
        CPPN(loser) = config.recFcn(CPPN(winner),CPPN(loser),config);
        CPPN(loser) = config.mutFcn(CPPN(loser),config);
        
        % Manipulate substrate
        substrate(loser) = config.recFcn(substrate(winner),substrate(loser),config_sub);
        substrate(loser) = config.mutFcn(substrate(loser),config_sub);
        
        %% Evaluate and update fitness
        [substrate(loser),~,CPPN(loser),~] = assessCPPNonSubstrate(substrate(loser),config_sub,CPPN(loser),config);
        
        %update errors
        store_error(test,gen,:) =  store_error(test,gen-1,:);
        store_error(test,gen,loser) = getError(config.error_to_check,substrate(loser));
        %genotype(loser).valError;
        best(gen)  = best(gen-1);
        best_indv(gen) = best_indv(gen-1);
        
        % print info
        if (mod(gen,config.gen_print) == 0)
            [best(gen),best_indv(gen)] = min(store_error(test,gen,:));
            fprintf('Gen %d, time taken: %.4f sec(s)\n  Winner: %.4f, Loser: %.4f, Best Error: %.4f \n',gen,toc/config.gen_print,getError(config.error_to_check,substrate(winner)),getError(config.error_to_check,substrate(loser)),getError('test',substrate(best_indv(gen))));
            tic;
            if best_indv(test,gen) ~= last_best
            % plot reservoir structure, task simulations etc.
            plotReservoirDetails(substrate,best_indv(test,:),gen,loser,config_sub); 
            plotReservoirDetails(CPPN,best_indv(test,:),gen,loser,config); 
            end
            last_best = best_indv(gen);    
        end
        
        %get metric details
        if config.record_metrics
            parfor pop_indx = 1:config.pop_size
                metrics(pop_indx,:) = getMetrics(genotype(pop_indx),config);
            end
        end
        
    end
end

function saveData(population,store_error,tests,config)
        config.figure_array =[];
        save(strcat('CPPN_substrate_',config.res_type,'_run',num2str(tests),'_gens',num2str(config.total_gens),'_',num2str(sum(config.num_reservoirs)),'Nres.mat'),...
            'population','store_error','config','-v7.3');
end