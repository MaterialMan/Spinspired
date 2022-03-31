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
config.pop_size = 30;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)
config.total_gens = 500;                    % number of generations to evolve

%% sweep film sizes
depth = round(linspace(0,5,config.pop_size)*10)/10;

config.test = 1;

% create population of reservoirs
population = config.createFcn(config);
default_pop = population(5);

ppm = ParforProgMon('Initial population: ', config.pop_size);
parfor pop_indx = 1:config.pop_size
    warning('off','all')
    tic;
    population(pop_indx) = default_pop;
    population(pop_indx).core_indx = pop_indx;
    population(pop_indx).damping = 1;
    population(pop_indx).system_size(3) = depth(pop_indx);
    population(pop_indx).behaviours = getMetrics(population(pop_indx),config);
    time_taken(pop_indx) = toc;
    ppm.increment();
end

% thickness vs time
figure
plot(depth,time_taken)
xlabel('Thickness (nm)')
ylabel('Seconds')

% metrics vs thickness
figure
plot(depth',reshape([population.behaviours],3,config.pop_size)')
%plot(repmat(time_taken,3,1)',reshape([population.behaviours],3,20)')
legend('KR','GR','MC')
xlabel('Thickness')
ylabel('Metric Value')