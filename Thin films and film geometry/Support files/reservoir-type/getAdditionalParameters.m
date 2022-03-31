%% Define additional params for particular reservoirs and tasks
% overflow for params that can be changed. This is called by all main scripts
function [config] = getAdditionalParameters(config)

%% Set Default parameters
config.mutate_type = 'gaussian';            %options: 'gaussian', 'uniform'. Type of distribution new weight is chosen from.
config.num_reservoirs = length(config.num_nodes);% num of subreservoirs. Default ESN should be 1.
config.leak_on = 1;                           % add leak states

% define connectivities
config.add_input_states = 1;                  %add input to states
config.sparse_input_weights = 1;              % use sparse inputs
config.sparsity = 0.1;                          % sparsity of input weights

config.internal_sparsity = 0.1;            % sparsity of internal reservoir  weights
config.input_weight_initialisation = 'norm';     % e.g.,  'norm', 'uniform', 'orth', etc. must be same length as number of subreservoirs
config.connecting_sparsity = 0;
config.internal_weight_initialisation = 'norm';  % e.g.,  'norm', 'uniform', 'orth', etc.  must be same length as number of subreservoirs

config.evolve_output_weights = 0;             % evolve rather than train
config.output_weight_initialisation = 'norm';  % e.g.,  'norm', 'uniform', 'orth', etc.  must be same length as number of subreservoirs
config.output_connectivity = 0.1;
config.output_weight_scaler = 1;              % defines maximum/minimum weight value when evolving output weights

config.evolve_feedback_weights = 0;             % find suitable feedback weights
config.feedback_weight_initialisation = 'norm';
config.feedback_connectivity = 0.1;
config.teacher_forcing = 0;                 % train output using target signal then transition into "generative" mode
config.feedback_scaling = 1;
config.noise_ratio = 0;                     % noise added in feedback training

% node functionality
config.activ_list = {@linearNode,@tanh};     % what activations are in use
config.multi_activ = length(config.activ_list) > 1;                      % use different activation funcs
config.training_type = 'Ridge';              % blank is psuedoinverse. Other options: Ridge, Bias,RLS
config.undirected = 0;                       % by default all networks are directed
config.undirected_ensemble = 0;              % by default all inter-network weights are directed

% default reservoir input scale
config.scaler = 1;                          % this may need to change for different reservoir systems that don't fit to the typical neuron range, e.g. [-1 1]
config.discrete = 0;                        % select '1' for binary input for discrete systems            
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems

% input mechanism to use
config.input_mechanism = 'continuous';         % if spiking, data is converted into spike trains. Otherwise uses normal analogue signals.

% default preprocessing performed on input data
config.preprocess = 'scaling';
config.preprocess_shift = [-1 1]; % range for data

% simulation details
config.run_sim = 0;
config.film = 0;                            %record simulation

% iir filter params (not currently in use)
config.iir_filter_on = 0;                       % add iir filters to states **not in use yet
config.iir_filter_order = 2 + 1;

%% Change/add parameters depending on reservoir type
% This section is for additional parameters needed for different reservoir
% types. If something is not here, it is most likely in the
% reservoir-specific @createFcn
if iscell(config.res_type)
    res_type = 'Heterotic';
else
    res_type = config.res_type;
    %config.figure_array = [figure figure];
end

switch(res_type)
    
    case {'RoR','RoRmin'}
        
        config.noise_level = 10e-6;
        config.mut_rate_connecting = 0.001;
        config.prune_rate = 0.00;
        
        config.RoR_structure = 0;
        if config.RoR_structure
            config.graph_type = {'Ring'};
            % Define substrate. Add graph type to cell array for multi-reservoirs
            % Examples: 'Hypercube','Cube'
            % 'Torus','L-shape','Bucky','Barbell','Ring',
            % 'basicLattice','partialLattice','fullLattice','basicCube','partialCube','fullCube',ensembleLattice,ensembleCube,ensembleShape
            config.self_loop = [0];               % give node a loop to self. Must be defined as array.
            % node details and connectivity
            [config,config.num_nodes] = getShape(config);
        end
        
    case 'ELM'
        config.leak_on = 0;                           % add leak states
        
    case 'BZ'
        %config.plot_BZ =0;
        config.fft = 0;
        config.evolve_output_weights = 0;             % evolve rather than train
        config.plot_states = 0;
        
        %input params
        config.sparse_input_weights = 1;
        config.input_widths = 1;
        config.bias_node = 0;
        config.max_time_period = 10;
        config.input_weight_initialisation = 'uniform';
        config.preprocess = 'scaling';
        config.preprocess_shift = 'zero to one';
        
    case 'Graph'
        
        
        config.graph_type= {'fullLattice'};            % Define substrate. Add graph type to cell array for multi-reservoirs
        % Examples: 'Hypercube','Cube'
        % 'Torus','L-shape','Bucky','Barbell','Ring',
        % 'basicLattice','partialLattice','fullLattice','basicCube','partialCube','fullCube',ensembleLattice,ensembleCube,ensembleShape
        config.self_loop = [0];               % give node a loop to self. Must be defined as array.
        
        % error checks
        if length(config.graph_type) ~= length(config.num_nodes) || length(config.self_loop) ~= length(config.num_nodes)
            error('Number of graph types does not match number of reservoirs. Add more.')
        end
        if contains(config.graph_type,'Lattice')
            if (sqrt(config.num_nodes) ~= round(sqrt(config.num_nodes)))
                error('\n Number of nodes needs to be a square number. \n')
            else
                config.num_nodes = sqrt(config.num_nodes);
            end
        end
        
        % node details and connectivity
        config.SW_type = 'topology';                                % maintain the base topology through evolution
        config.ensemble_graph = 0;                                  % no connections between mutli-graph reservoirs
        [config,config.num_nodes] = getShape(config);               % call function to make graph.
        
        % pick node functions
        config.multi_activ = 0;                      % use different activation funcs
        config.activ_list = {@tanh};     % what activations are in use
        
        
    case 'DNA'
        config.tau = 20;                         % settling time
        config.step_size = 1;                   % step size for ODE solver
        config.concat_states = 0;                % use all states
        
    case 'RBN'
        
        config.k = 2; % number of inputs
        config.mono_rule = 1; % use one rule for every cell/reservoir
        config.rule_list = {@evolveCRBN}; %list of evaluation types: {'CRBN','ARBN','DARBN','GARBN','DGARBN'};
        config.leak_on = 0;
        config.discrete = 1;
        config.time_period = 1;
        
    case 'elementary_CA'
        % update type
        config.k = 3;
        config.mono_rule = 1;               %stick to rule rule set, individual cells cannot have different rules
        config.rule_list = {@evolveCRBN}; %list of evaluation types: {'CRBN','ARBN','DARBN','GARBN','DGARBN'};
        config.leak_on = 0;
        config.discrete = 1;
        config.torus_rings = 1;
        config.rule_type = 0;
        
    case '2D_CA'
        % update type
        config.sparse_input_weights = 1;
        config.mono_rule = 1;               %stick to rule rule set, individual cells cannot have different rules
        config.rule_list = {@evolveCRBN}; %list of evaluation types: {'CRBN','ARBN','DARBN','GARBN','DGARBN'};
        config.leak_on = 0;
        config.rule_type = 'Moores';
        config.discrete = 1;
        config.torus_rings = 1;
        
    case 'DL'
        %     config.DLtype = 'mackey_glass2';%'ELM';%'virtualNodes';
        %     %config.tau = 100;
        config.preprocess = 0;
        config.tau = config.num_nodes.*0.2; % keep 0.2 separation at all sizes
        config.binary_weights = 0;
        
    case 'CNT'
        
        config.volt_range = 5;
        config.num_input_electrodes = 64;
        config.num_output_electrodes = 32;
        config.discrete = 0;
        
    case 'Wave'
        
        config.leak_on = 1;                           % add leak states
        config.add_input_states = 1;                  %add input to states
        config.input_widths = 0;
        config.sim_speed = 1; % xfactor
        config.time_step = 0.05;
        config.bias_node = 0;
        config.max_time_period = 20;
        config.max_wave_speed = 12;
        % [1 0 0] = fix: All boundary points have a constant value of 1
        % [0 1 0] = cont; Eliminate the wave and bring elements to their steady state.
        % [0 0 1] = connect; Water flows across the edges and comes back from the opposite side
        config.boundary_conditions = {[1 0 0],[0 1 0],[0 0 1]}; %add more as cell, if needed
        config.wave_system = 'ensemble';% type of multi-res system. Choices "ensemble", "fully-connected", "pipeline"
        
        if length(config.num_nodes) == 1
            config.wave_system = 'ensemble'; % no effect on system if one material
        end
        
        
    case 'GOL'
        config.leak_on = 1;                           % add leak states
        config.add_input_states = 1;                  %add input to states
        config.sparse_input_weights = 1;
        config.discrete = 0;
        
    case 'SW'
        
        config.graph_type= {'Ring'};            % Define substrate. Add graph type to cell array for multi-reservoirs
        % Examples: 'Hypercube','Cube'
        % 'Torus','L-shape','Bucky','Barbell','Ring',
        % 'basicLattice','partialLattice','fullLattice','basicCube','partialCube','fullCube',ensembleLattice,ensembleCube,ensembleShape
        config.self_loop = [1];               % give node a loop to self. Must be defined as array.
        
        if length(config.graph_type) ~= length(config.num_nodes) && length(config.self_loop) ~= length(config.num_nodes)
            error('Number of graph types does not match number of reservoirs. Add more in getDataSetInfo.m')
        end
        
        % node details and connectivity
        config.ensemble_graph = 0;                                  % no connections between mutli-graph reservoirs
        config.P_rc = 0.1;                                         % percentage of random connections or connection probability. Used for Small World Networks
        [config,config.num_nodes] = getShape(config);               % call function to make graph.
        
        config.SW_type = 'topology_plus_weights';                    %options: 'topology_plus_weights', 'watts_strogartz'
        
    case 'Pipeline'
        config.add_input_states = 1;
        
    case 'Oregonator'
        
        config.add_input_states = 1;
        config.figure_array = [config.figure_array figure];
        config.sparse_input_weights = 1;
        config.sparsity(1) = 0.2;
        config.input_weight_initialisation = 'uniform';
        config.leak_on = 1;
        config.bias_node = 0;
        
        config.stride      = 100;
        config.displayLive = true;
        config.displayZoom = 4.0;
        
        config.epsilon     = 0.0243; %0.0243
        config.grid_2      = 0.0625; % 0.25^2
        config.diff_coeff  = 0.45;
        config.f           = 1.4;
        config.q           = 0.002; %0.002
        config.speed       = 0.0005;
        config.phi_active  = 0.054;    % normal, excitable
        config.phi_passive = 0.0975;   % passive, excitable
        
        config.u_idx = 1;
        config.v_idx = 2;
        config.p_idx = 3;
        config.b_idx = 4;
        
        config.crossover_dist = 1;
        
        config.vesicle_radius = 5;
        %config.min_time_period = 100;
        config.max_time_period = config.stride;
        %config.input_length = 10;
        config.plot_states = 0;
        
        % must be non-negative
        config.preprocess = 'scaling';
        config.preprocess_shift = 'zero to one';
        
    case 'Heterotic'
        % reservoir params
        
        % multi-reservoir type
        config.architecture = 'pipeline_IA'; % can be 'ensemble','pipeline', or 'RoR'
        config.plot_states = 0;
        
        % error checks
        if length(config.res_type) ~= length(config.num_nodes)
            error('Number of graph types does not match number of reservoirs. Add more.')
        end
        
case {'MM','multiMM'}
        % default MM data scaling
        config.preprocess = 'rescale_diff';
        config.preprocess_shift = [0 1]; % range for data
        
        % reservoir params
        config.leak_on = 1;                         % add a leak rate filter
        config.add_input_states = 1;                % add input to states
        config.single_input = 1;        
        config.bias = 1;                            % whether to apply a bias to cells as an additional input; bias = value given, weights then alter this for different cells
        config.sparse_input_weights = 1;            % use a sparse input encoding
        config.sparsity = 0.1;                    % 0 to 1 sparsity of input weights
        config.input_widths = 0;                    % inputs can apply to a single cell (when '0') or multiple cells with a disspating radius of 'r' (r>0)
        config.input_scaler = 2;                    % input weight multiplier, e.g.x2
        config.input_weight_initialisation = 'norm';     % e.g.,  'norm', 'uniform', 'orth', etc. Check createMM.m for options
        
        % system settings
        config.material_type = {'toy'};   % options: 'toy', 'multilayer','core_shell', 'random_alloy', '' (if specific config)
        config.crystal_structure = {'sc'};            % typical crystal structures: 'sc', 'fcc', 'bcc' | 'sc' significantly faster
        config.unit_cell_size = [2.507];               % depends on crystal structure; typical value 3.47 Armstrongs fo 'sc'
        config.unit_cell_units = {'!A'};              % range = 0.1 � to 10 � m
        config.macro_cell_size = [5];                % size of macro cell; an averaging cell over all spins inside
        config.macro_cell_units = {'!nm'};            % units for macro cell size
        
        config.system_size_z = [0.1];               % thickness of film/system; currently x and y are determined by node size and is always a square. *Macro size needs to be twice the thickness Z(?)
        config.size_units = {'!nm','!nm','!nm'};     % units for x,y,z size
        
        % additional properties: random alloys, core shells, periodic
        % boundaries, user-specific strucutres etc.
        config.particle_size = [];                  % only in use when used with: shapes, core shell ; must be less than system size!
        config.periodic_boundary = [0,0,0];         % vector represents x,y,z; '1' means there is a periodic boundary
        config.material_shape = {'film'};             % type shape to cut out of film; check shape is possible,e.g. film is default
                
        config.evolve_geometry = 1;                    % manipulate geomtry
        config.evolve_poly = 0; % otherwise evolve a rectangle
        config.poly_num = 4; 
        config.geometry_file =  'custom.geo';                    %add specific geometry file: '' is blank (when not in use), use 'custom.geo' to evolve a polygon
        config.lb = 0.1;                                %lower bound on dimension

        %defaults
        for i = 1:length(config.num_nodes)
            config.random_alloy(i) = false;
            config.core_shell(i) = false;
            config.user_structure_file{i} = '';            % if wanting to load a specific strucutre, write file name
            config.evolve_material_density(i) = false;     % whether density can be changed: 'static' or 'dynamic'
            
            switch(config.material_type{1})
                case 'toy'
                    config.unit_cell_size = 3.47;
                    config.crystal_structure = {'sc'}; 
                    config.unit_cell_units = {'!A'};
                case 'random_alloy'
                    % specific to fcc structure
                    config.unit_cell_size = 3.524;
                    config.unit_cell_units = {'!A'};
                    config.random_alloy(i) = true;                % set to random alloy
                case 'core_shell'
                    config.core_shell(i) = true;                  % apply a core shell configuration
                case 'multilayer'
                    
                case 'user_defined'
                    config.user_structure_file = {''};
            end
        end
        config.shell_sizes = [];                    % defines the radial extent of a material as a fraction of the particle radius, e.g., [shell, core], range = 0 to 1
        %config.material_heights = [];               % type end heights of each material; if more than one, will start from zero to height(1), then height(1) to height(2) and so on.
        %config.evolve_material_density = [false];     % whether density can be changed: 'static' or 'dynamic'
        config.material_density = [1];              % removes random atoms to equal density; 0 to 1 for each material, e.g., [0.5, 1]
        
        % material properties
        config.temperature_parameter = [0,0];           % positive integer OR 'dynamic'
        config.temperature_rescaling_exponent = [2.369 2.369]; % for more acurate temperature measures. [Co=2.369,FE=2.876,Ni=2.322]
        config.temperature_rescaling_curie_temperature = [1395 1395]; % [Co=1395,FE=1049,Ni=635]
        config.damping_parameter = [0.001, 1];             % 0 to 10 OR 'dynamic' | typical value 0.1
        config.anisotropy_parameter = [6.69e-24, 6.69e-24];   % 1e-26 to 1e-22 OR 'dynamic' | typical value 1e-24
        config.exchange_parameter = [6.064e-21, 6.064e-21];    % 1e-21 to 10e-21 OR 'dynamic' | typical value 5e-21
        config.magmoment_parameter = [1.72, 1.72];            % 1 (<1muB can have intergration problems) to 10 OR 'dynamic' | typical value 1.4
        config.applied_field_strength = [0,0];          % how many tesla (T)
        config.initial_spin_direction = {'1,0,0'};      % assign initial spin direction of material as a string. Add more cells for more materials.
        
        config.interfacial_exchange = [-25e-21, 25e-21];
        
        %simulation params
        config.time_step = 100;                    % simulation/itegrator time-step
        config.time_units = '!fs';                  % must have '!' before unit
        config.time_steps_increment = [100 100];    % time step to apply input; e.g. 100 or 1000
        config.read_mag_direction = {'z'};  % list of directions to read; can be 1, 2  or all
        config.applied_field_unit_vector = {'0,0,1'}; % where the applied field will be directed; x,y,z
        
        % plot output
        config.plt_system = 0;                      % to create plot files from vampire simulation
        config.plot_rate =100;                        % rate to build plot files
        config.plot_states = 0;                     % plot every state in matlab figure; for debugging
        
        % multi-reservoir type
        config.architecture = '';
             
        % calculate number of reservoirs in total if multilayered system (signified by cell type)
        if iscell(config.num_nodes)
            config.num_layers = length(config.num_nodes);
            for i = 1:config.num_layers 
                len(i) = length(config.num_nodes{i});
            end
            config.dummy_node_list = zeros(config.num_layers,max(len));
            
            for i = 1:config.num_layers  % cycle through layers
                config.dummy_node_list(i,1:length(config.num_nodes{i})) = config.num_nodes{i};
                config.num_res_in_layer(i) = length(config.num_nodes{i});
                config.total_units_per_layer(i) = sum(config.num_nodes{i});
            end
            config.total_units = sum(config.total_units_per_layer);
        end
        
       % params for multilayered system
        config.add_pipeline_input = 0;
        %config.num_layers = size(config.num_nodes,1); % calculates layers from cell structure in num_nodes
        config.sparsity = 0.1; % input sparsity
        config.input_weight_initialisation = 'norm';     % e.g.,  'norm', 'uniform', 'orth', etc. must be same length as number of subreservoirs
        config.connecting_sparsity = 0.001; % connecting sparsity
        config.internal_weight_initialisation = 'norm';  % e.g.,  'norm', 'uniform', 'orth', etc.  must be same length as number of subreservoirs
        
    otherwise
        
end

%% Task parameters - now apply task-specific parameters
% If a task requires additional parameters, or resetting from defaults, add
% here.
switch(config.dataset)
    
    case 'cifar10'
        config.leak_on = 0;                          % add leak states
        %config.add_input_states = 0;
        
    case 'autoencoder'
        config.leak_on = 0;                          % add leak states
        config.add_input_states = 0;
        config.figure_array = [config.figure_array figure];
        config.sparse_input_weights = 1;
        %config.evolve_output_weights = 1;               % must evolve outputs as its an unsupervised problem
        %config.error_to_check = 'train&val';
        config.preprocess = 'clip';
        config.preprocess_shift = '';
        
    case 'pole_balance'
        config.time_steps = 1000;                        % length of task simulation
        config.simple_task = 2;                         % tasks: 1) balancing pole from up-right position , 2) swinging pole from downward position , 3)
        config.pole_tests = 3;                          % how many tasks to average over
        config.velocity = 1;                            % add velocity as an input to task (usually easier)
        config.run_sim = 0;                             % whether to run simulation as it is calculated
        config.testFcn = @poleBalance;
        config.evolve_output_weights = 1;               % must evolve outputs as its an unsupervised problem
        config.evolve_feedback_weights = 0;
        config.leak_on = 1;
        config.add_input_states = 0;                    % add input to states
        config.error_to_check = 'train';
        
    case 'robot'
        % type of task
        config.robot_behaviour = 'explore_maze';    %select behaviour/file to simulate
        config.time_steps = 500;                    % sim time
        %sensors
        config.sensor_range = 0.5;                 % range of lidar
        config.evolve_sensor_range = 1;             % use leakRate parameter as proxy for sensor range (evolvable)
        config.sensor_radius = 2*pi;
        % sim parameters
        config.run_sim = 0;                          % whether to run/watch sim
        config.robot_tests = 3;                     % how many tests to conduct: to provide avg fitness
        config.show_robot_tests = 1; % how many tests to watch/check visually
        config.sim_speed = 30;                       % speed of sim result/visualisation. e.g. if =2, 2x speed
        config.testFcn = @robot;                    % assess fcn for robot tasks
        config.evolve_output_weights = 1;             % must be on; unsupervised/reinforcement problem
        
        %environment
        config.bounds_x = 10;                        % scaler for extending bounds of environment
        config.bounds_y = 10;
        config.num_obstacles = 10;                   % number of obstacles to place in environment
        config.num_target_points = 1000;            % grid of target points used for fitness calculation
        config.maze_size = 5;                       % if maze, the size and complexity of maze
        % Go to selectDataset.m to change num_sensors
        config.error_to_check = 'train';
        
    case 'attractor'
        
        config.error_to_check = 'train&test';
        config.attractor_type = 'mackey_glass';
        config.preprocess = '';
        config.preprocess_shift = [0 1];

        config.leak_on = 0;                          % remove leak states
        config.add_input_states = 0;                 % remove input states
        config.sparse_input_weights = 1;
        
        config.teacher_forcing = 1;                 % train output using target signal then transition into "generative" mode; evolving output weights will be worse
        config.noise_ratio = 10e-5;                     % noise added in feedback training

        config.evolve_output_weights = 1;             % evolve rather than train
        config.output_weight_initialisation = 'uniform';  % e.g.,  'norm', 'uniform', 'orth', etc.  must be same length as number of subreservoirs
        config.output_connectivity = 0.1;
        config.output_weight_scaler = 1;              % defines maximum/minimum weight value when evolving output weights
        
        config.evolve_feedback_weights = 1;             % find suitable feedback weights - currently doesn't work with teacher forcing
        config.feedback_weight_initialisation = 'uniform';
        config.feedback_connectivity = 0.1;
        config.feedback_scaling = 1;

    case {'MSO1','MSO2','MSO3','MSO4','MSO5','MSO6','MSO7','MSO8','MSO9','MSO10','MSO11','MSO12'}  %MSO'
        
        %config.preprocess = 0;
        
        
    case {'image_painting','CPPN'}
        %config.leak_on = 0;                          % add leak states
        config.add_input_states = 0;
        
        % dummy figure
        %config.figure_array = [];
        
    case{'test_pulse'}
        config.preprocess = 0;
        config.preprocess_shift = [0 1];
        
    otherwise
        
end


