
clear %vars -except population
close all

rng(1,'twister');

%% Setup
config.parallel = 1;                        % use parallel toolbox

%start paralllel pool if empty
if isempty(gcp) && config.parallel
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

res_list = [25,100,225,400];

for res_size = 1:4
    % type of network to evolve
    config.res_type = 'MM';            % state type of reservoir(s) to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network with multiple functions), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more. Place reservoirs in cell ({}) for heterotic systems.
    config.num_nodes = [res_list(res_size)];                   % num of nodes in each sub-reservoir, e.g. if config.num_nodes = [10,5,15], there would be 3 sub-reservoirs with 10, 5 and 15 nodes each.
    config = selectReservoirType(config);         % collect function pointers for the selected reservoir type
    
    %% Evolutionary parameters
    config.pop_size = 50;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)
    config.error_to_check = 'train&val&test';
    
    %% Task parameters
    config.dataset = 'laser';          % Task to evolve for
    config.figure_array = [figure figure];
    
    % get any additional params. This might include:
    % details on reservoir structure, extra task variables, etc.
    config = getAdditionalParameters(config);
    
    % get dataset information
    config = selectDataset(config);
    
    %% general params
    % type of metrics to apply; if necessary
    config.metrics = {'KR','GR','linearMC'};          % list metrics to apply in cell array: see getVirtualMetrics.m for types of metrics available
    config.record_metrics = 0;                  % save metrics
    
    %set random seed
    rng(1,'twister');
    config.test = 1;
    
    % create initial population
    population = config.createFcn(config);
    
    is_range = [0.8 0.95];
    d_range =[0.175 0.3];
    lr_range =[0.7 0.9];
    
    for i = 1:config.pop_size
        population(i).input_scaling = is_range(1) + (is_range(2)-is_range(1))*rand;
        population(i).damping = d_range(1) + (d_range(2)-d_range(1))*rand;
        population(i).leak_rate = lr_range(1) + (lr_range(2)-lr_range(1))*rand;
    end
    
    %Assess population
    ppm = ParforProgMon('Initial population: ', config.pop_size);
    parfor pop_indx = 1:config.pop_size
        warning('off','all')
        population(pop_indx) = config.testFcn(population(pop_indx),config);
        fprintf('\n i = %d, error = %.4f\n',pop_indx,getError(config.error_to_check,population(pop_indx)));
        ppm.increment();
    end
    
    store_pop(res_size,:) = population;
end