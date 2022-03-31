
%% Evolve substrate for a specific task
% This script can be used to evolve any reservoir directly to a task. It
% uses the steady-state Microbial Genetic Algorithm to evolve the best
% solution.

% Author: M. Dale
% Date: 03/07/19
function population = STMFcn(ID,res_type,node_size,task,pop_size,num_gens,hetero_exchange_interaction,array_type)

%delete(gcp('nocreate'))
%parpool('local');

addpath(genpath('/mnt/lustre/users/md596/Magnets/multilayer'));

rng(ID,'twister');

%% Setup
config.test = ID;
config.parallel = 0;                        % use parallel toolbox

% type of network to evolve
config.res_type = res_type;             % state type of reservoir(s) to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network with multiple functions), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = node_size;                   % num of nodes in each sub-reservoir, e.g. if config.num_nodes = [10,5,15], there would be 3 sub-reservoirs with 10, 5 and 15 nodes each.
config = selectReservoirType(config);         % collect function pointers for the selected reservoir type

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = pop_size;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)
config.total_gens = num_gens;                    % number of generations to evolve
config.mut_rate = 0.01;                       % mutation rate
config.deme_percent = 0.2;                   % speciation percentage; determines interbreeding distance on a ring.
config.deme = round(config.pop_size*config.deme_percent);
config.rec_rate = 0.5;                       % recombination rate
config.error_to_check = 'train&val&test';    % printed error includes all three errors

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = task;          % Task to evolve for

config.hetero_exchange_interaction = hetero_exchange_interaction;

% get any additional params. This might include:
% details on reservoir structure, extra task variables, etc.
config = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% general params
config.gen_print = 10;                       % after 'gen_print' generations print task performance and show any plots
config.start_time = datestr(now, 'HH:MM:SS');
config.save_gen = config.gen_print;                       % save data at generation = save_gen
config.figure_array = [figure figure];

% Only necessary if wanting to parallelise the microGA algorithm
config.multi_offspring = 0;                 % multiple tournament selection and offspring in one cycle
config.num_sync_offspring = 0;%config.deme;    % length of cycle/synchronisation step

% type of metrics to apply; if necessary
config.metrics = {'linearMC'};          % list metrics to apply in cell array: see getVirtualMetrics.m for types of metrics available
config.record_metrics = 1;                  % save metrics

% additional pruning, if required
config.prune = 0;                           % after best individual is found prune the network to make more efficient
config.tolerance = [0 0 0];
config.prune_iterations = 500;

%% Run experiments
for test = 1:config.num_tests
    
    clearvars -except config test best best_indv store_error population ID node_size task array_type
    
    warning('off','all')
    fprintf('\n Test: %d  ',test);
    fprintf('Processing genotype......... %s, task is %s \n',datestr(now, 'HH:MM:SS'),task)
    tic
    
    % tracker for sim plot
    last_best = 1;
    
    %set random seed
    rng(ID + test,'twister');
    
    % change array struct
     config.graph_type = {array_type};
     [~,~,config.G{1}] = getShape(config,config.num_nodes{1},'diagraph');

    % create initial population
    population = config.createFcn(config);
    

    % apply metrics to final population
    if config.record_metrics
        parfor pop_indx = 1:config.pop_size
            population(pop_indx).core_indx = config.test + pop_indx;
            input_weights = zeros(2,population(pop_indx).layer(1).nodes);
            input_weights(:,1) = 1;
            population(pop_indx).layer(1).input_weights{1}=input_weights;
            
            metrics(pop_indx,:) = getMetrics(population(pop_indx),config);
            population(pop_indx).metrics = metrics(pop_indx,:);
        end
    end
    
end
end

