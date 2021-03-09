function population = createRestrictedRoR(config)

%% Reservoir Parameters
for pop_indx = 1:config.pop_size
    
    % add performance records
    population(pop_indx).train_error = 1;
    population(pop_indx).val_error = 1;
    population(pop_indx).test_error = 1;
    
    % add single bias node
    population(pop_indx).bias_node = 1;
    
    % assign input/output count
    if isempty(config.train_input_sequence)
        population(pop_indx).n_input_units = 1;
        population(pop_indx).n_output_units = 1;
    else
        population(pop_indx).n_input_units = size(config.train_input_sequence,2);
        population(pop_indx).n_output_units = size(config.train_output_sequence,2);
    end
    
    
    %define num of units
    population(pop_indx).nodes = config.num_nodes;
    population(pop_indx).total_units = sum(config.num_nodes);
    
    % Scaling and leak rate
    population(pop_indx).init_input_scaling = 2*rand(1,length(config.num_nodes))-1; % can expand to have for each input and bias: population(pop_indx).n_input_units+1
    population(pop_indx).init_W_scaling = 2*rand(1,length(config.num_nodes))-1;
    population(pop_indx).init_leak_rate = rand(1,sum(config.num_nodes));
    
    %inputweights
    population(pop_indx).input_weights = getWeights(config.input_weight_initialisation...
        ,population(pop_indx).n_input_units+1, population(pop_indx).total_units,config.sparsity);
    
    % restrict number of inputs
    population(pop_indx).restricted_input_mask = zeros(population(pop_indx).n_input_units+1, population(pop_indx).total_units);
    input_pos = randperm(population(pop_indx).total_units*(population(pop_indx).n_input_units+1), config.restricted_num_inputs);%randi([1 population(pop_indx).nodes(i)*(population(pop_indx).n_input_units+1)], config.restricted_num_inputs,1);
    population(pop_indx).restricted_input_mask(input_pos) = 1;
    
    %assign different activations, if necessary
    if config.multi_activ
        population(pop_indx).activ_Fcn_indx = randi(length(config.activ_list),1,population(pop_indx).total_units);
    end
    
    population(pop_indx).last_state = zeros(1,population(pop_indx).total_units);
    
    %% weights and connectivity of all reservoirs
    % end of first subres
    end_pos = population(pop_indx).nodes(1);
    
    % intialise W and W_scale
    population(pop_indx).W_scaling(1:end_pos,1:end_pos) = population(pop_indx).init_W_scaling(1);
    population(pop_indx).input_scaling(:,1:end_pos) = population(pop_indx).init_input_scaling(:,1);
    population(pop_indx).leak_rate(1:end_pos) = population(pop_indx).init_leak_rate(1);
    
    % assign first internal weights
    population(pop_indx).W = getWeights(config.internal_weight_initialisation,...
        population(pop_indx).nodes(1),population(pop_indx).nodes(1),config.internal_sparsity);
    
    for i = 2:length(population(pop_indx).nodes)
        
        start_pos = end_pos+1;
        end_pos = start_pos + population(pop_indx).nodes(i)-1;
        
        population(pop_indx).input_scaling(:,start_pos:end_pos) = population(pop_indx).init_input_scaling(:,i);
        population(pop_indx).W_scaling(start_pos:end_pos,start_pos:end_pos) = population(pop_indx).init_W_scaling(i);
        
        population(pop_indx).leak_rate(start_pos:end_pos) = population(pop_indx).init_leak_rate(i);
        
        % assign new subres weights
        population(pop_indx).W(start_pos:end_pos,start_pos:end_pos) = getWeights(config.internal_weight_initialisation,...
            population(pop_indx).nodes(i),population(pop_indx).nodes(i),config.internal_sparsity);
        
    end
    
    % mask to indicate separation between reservoirs
    population(pop_indx).test_mask = (population(pop_indx).W_scaling == 0);
    
    % Architecture switch - if defined architecture is given for RoR
    % system, e.g. ring or lattice of reservoirs
    if config.RoR_structure
        % find indices for graph weights
        graph_indx = logical(full(adjacency(config.G{1})));
        % assign connectivity
        population(pop_indx).W_switch = graph_indx;
    else
        population(pop_indx).W_switch = round(rand(config.num_reservoirs));
    end
    
    
    % add rand output weights
    if config.add_input_states
        output_units = population(pop_indx).total_units + population(pop_indx).n_input_units;
    else
        output_units = population(pop_indx).total_units;
    end
    
    % created restricted output mask
    population(pop_indx).restricted_output_mask = zeros(population(pop_indx).total_units, 1);
    output_pos = randperm(population(pop_indx).total_units,config.restricted_num_outputs);%randi([1 population(pop_indx).nodes(i)],config.restricted_num_outputs, 1);
    population(pop_indx).restricted_output_mask(output_pos) = 1;
    
    %apply initialisation technique
    population(pop_indx).output_weights = getWeights(config.output_weight_initialisation...
        ,output_units, population(pop_indx).n_output_units,config.output_connectivity);
    
    % add rand feedback weights
    if config.evolve_feedback_weights
        population(pop_indx).feedback_scaling = 2*rand-1;
        population(pop_indx).feedback_weights = getWeights(config.feedback_weight_initialisation...
            ,population(pop_indx).total_units, population(pop_indx).n_output_units,config.output_connectivity);
    end
    
    % assign blank space for behaviours for CHARC characterisation
    population(pop_indx).behaviours = [];
end
end

function weights = getWeights(weight_initialisation,x_dim,y_dim,sparsity)

switch(weight_initialisation)
    case 'norm' % normal distribution
        weights = sprandn(x_dim, y_dim, sparsity);
    case 'uniform' % uniform dist between -1 and 1
        weights = sprand(x_dim, y_dim, sparsity);
        weights(weights ~= 0) = ...
            2*weights(weights ~= 0)  - 1;
    case 'orth'
        weights = orth(rand(x_dim, y_dim));
    case 'sparse_orth'
        weights = orth(full(sprand(x_dim, y_dim, sparsity)));
end
end