
function shape_error = shapeTest(ID,res_size,dataset,pop_size)

close all

addpath(genpath('/mnt/lustre/users/md596/Magnets/multilayer'));

%set random seed for experiments
rng(ID,'twister');

% type of network to evolve
config.res_type = 'multiMM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = {[res_size]};                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = dataset;          % Task to evolve for
config.figure_array = [figure figure figure];
config.preprocess ='';

config.metrics = {'KR','GR','linearMC'};

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = pop_size;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = ID;

%% shape parameters
%base_size = 64;
%config.total_units = base_size;

shape_list = {'parallelogram','trapezoid','rectangle','custom','square','custom','circle'};
name_list = {'0','1','2','3','4','5','inf'};
shape_num_sides = [4, 4, 4, 3, 4, 5, 40];
shape_in_pos = {46,[46,54],[48,53],[53,54],[45,46,55,56],51,[39,40,49,50]};

%scale_range = {[],[],[],[],[],[],[]};
%damping_range = {[],[],[],[],[],[],[]};

config.rotate_angle = 0;
config.ref_point = [0.5 0.5];

config.add_input_states = 1;

%cycle through shapes
for shape_indx = 1:length(shape_list)
    
    config.shape_list = shape_list{shape_indx}; 
    config.shape_num_sides = shape_num_sides(shape_indx);
    
    % define individual
    default_pop = config.createFcn(config);   
              
    parfor pop_indx = 1:length(default_pop)

        rng(pop_indx,'twister');
        
        % apply inputs
        default_pop(pop_indx).layer.input_weights{1} = zeros(size(default_pop(pop_indx).layer.input_weights{1}));
        
        % view input location
        s  = zeros(sqrt(size(default_pop(pop_indx).layer.input_weights{1},2)));
        s(default_pop(pop_indx).layer.inputs_mask) = 1:length(default_pop(pop_indx).layer.inputs_mask);
        s(default_pop(pop_indx).layer.inputs_mask(shape_in_pos{shape_indx})) = 2;
        imagesc(s)
        
        %apply input to centre of shape
        default_pop(pop_indx).layer.input_weights{1}(:,default_pop(pop_indx).layer.inputs_mask(shape_in_pos{shape_indx})) = 1;       
        
        % reset params
        default_pop(pop_indx).layer.leak_rate = 1;
        %default_pop(pop_indx).layer.input_scaling = scale_range{shape_indx}(1) + (scale_range{shape_indx}(2)-scale_range{shape_indx}(1))*rand;
        %default_pop(pop_indx).layer.damping = damping_range{shape_indx}(1) + (damping_range{shape_indx}(2)-damping_range{shape_indx}(1))*rand;

        % evaluate
        default_pop(pop_indx) = config.testFcn(default_pop(pop_indx),config);
        
        fprintf('\n shape=%s , i = %d, error = %.4f\n',config.shape_list,pop_indx,getError('test',default_pop(pop_indx)));
    
        shape_error(shape_indx,pop_indx) = default_pop(pop_indx).test_error;
        
    end
    
    shape_pop{shape_indx} = default_pop;
    
    save(strcat('shape_results_',dataset,'_size',num2str(res_size),'_pop',num2str(pop_size),'_run',num2str(ID),'.mat'),'shape_error','shape_pop')
end