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

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = 'step_response';          % Task to evolve for

config.figure_array = [figure figure];

config.metrics = {'linearMC'}; 

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = 5;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = 1;

% create population of reservoirs
population = config.createFcn(config);
default_pop = population(5);

%% damping effect
damping = linspace(0,1,20);

ppm = ParforProgMon('Initial population: ', length(damping));
parfor pop_indx = 1:length(damping)
    warning('off','all')

    population(pop_indx) = default_pop;
    population(pop_indx).input_scaling = 1;
    population(pop_indx).leak_rate = 1;
    population(pop_indx).input_weights{1} = zeros(100,2);
    population(pop_indx).input_weights{1}(50,1) = 1;
    
    population(pop_indx).core_indx = pop_indx;
    population(pop_indx).damping = damping(pop_indx);

    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
    MC(pop_indx) = getMetrics(population(pop_indx),config);
    
    ppm.increment();
end

% thickness vs time
colors = distinguishable_colors(length(damping));
figure
hold on
for i = 4:length(damping)
    p = plot(states{i}(:,201:300),'Color',colors(i,:));
    for n = 2:100
        set(get(get(p(n),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
hold off
xlabel('Time')
ylabel('Magnitude')
[~,hobj] = legend(string(round(damping(4:length(damping))*100)/100));
h1 = findobj(hobj,'type','line');
set(h1,'LineWidth',4);

figure
bar(MC)
xticks(1:20)
xticklabels(round(damping*100)/100)
xtickangle(45)
xlabel('Damping')
ylabel('MC')