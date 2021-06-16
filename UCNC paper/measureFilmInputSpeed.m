clear
close all

%set random seed for experiments
rng(1,'twister');

%% Setup
config.parallel = 0;                        % use parallel toolbox

%start paralllel pool if empty
if isempty(gcp) && config.parallel %
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

% type of network to evolve
config.res_type = 'MM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = [16*16];                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = 'test_pulse';          % Task to evolve for
config.figure_array = [figure figure];
config.preprocess ='';

config.metrics = {'linearMC'}; 

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = 1;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = 1;

% create population of reservoirs
population = config.createFcn(config);
default_pop = population(1);


%% damping effect
damping = linspace(0,1,15);

ppm = ParforProgMon('Initial population: ', length(damping));
parfor pop_indx = 1:length(damping)
    warning('off','all')

    population(pop_indx) = default_pop;
    population(pop_indx).input_scaling = 1;
    population(pop_indx).leak_rate = 1;
    population(pop_indx).input_weights{1} = zeros(config.num_nodes,2);
    population(pop_indx).input_weights{1}(config.num_nodes/2,1) = 1;
    
    population(pop_indx).core_indx = pop_indx;
    population(pop_indx).damping = damping(pop_indx);

    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
    %MC(pop_indx) = getMetrics(population(pop_indx),config);
    
    ppm.increment();
end

figure
n = 5;
for t = 1:size(states{n},1)
   imagesc(reshape(states{n}(t,1:end-1),sqrt(config.num_nodes),sqrt(config.num_nodes)))
   caxis([min(min(states{n})) max(max(states{n}))])
   colormap(bluewhitered)
   drawnow
   pause(0.1)
end
figure
%plot(states{n}(:,[50 41 1 99]))
plot(states{n}(:,config.num_nodes/2 - sqrt(config.num_nodes)+1))

for i = 1:length(damping)
[~,t_loc] = findpeaks(states{i}(:,config.num_nodes/2 - sqrt(config.num_nodes)+1), 'NPeaks',1,'MinPeakProminence',max(states{i}(:,config.num_nodes/2 - sqrt(config.num_nodes)+1))/2)
loc(i) = t_loc(1);
end

% calculate speed
freq = 100e9;
num_samples_to_peak = loc-10;

time = (1/freq).*num_samples_to_peak;
distance = (population(1).system_size(1)+1 - population(1).macro_cell_size*2 )*1e-9;

speed = distance./time;
figure
plot(damping,speed)

figure
scatter(speed(damping > 0.1),MC(damping > 0.1),50,damping(damping > 0.1),'filled')
ylabel('MC')
xlabel('Speed (m/s)')
title(num2str(config.num_nodes))
colorbar