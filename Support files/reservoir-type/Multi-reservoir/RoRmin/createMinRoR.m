function population = createMinRoR(config)


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
    population(pop_indx).input_weights = getInternalWeights(config.input_weight_initialisation...
        ,population(pop_indx).n_input_units+1, population(pop_indx).total_units,config.sparsity);
    
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
    
    if config.mulit_leak_rate 
        population(pop_indx).leak_rate(1:end_pos) = rand(1,sum(config.num_nodes)); % assign indvidual leak rates per node
    else
        population(pop_indx).leak_rate(1:end_pos) = population(pop_indx).init_leak_rate(1);
    end
    
    % assign first internal weights
    population(pop_indx).W = getInternalWeights(config.internal_weight_initialisation,...
        population(pop_indx).nodes(1),population(pop_indx).nodes(1),config.internal_sparsity);
    
    for i = 2:length(population(pop_indx).nodes)
        
        start_pos = end_pos+1;
        end_pos = start_pos + population(pop_indx).nodes(i)-1;
        
        population(pop_indx).input_scaling(:,start_pos:end_pos) = population(pop_indx).init_input_scaling(:,i);
        population(pop_indx).W_scaling(start_pos:end_pos,start_pos:end_pos) = population(pop_indx).init_W_scaling(i);
        
        if ~config.mulit_leak_rate
            population(pop_indx).leak_rate(start_pos:end_pos) = population(pop_indx).init_leak_rate(i);
        end
        
        % assign new subres weights
        population(pop_indx).W(start_pos:end_pos,start_pos:end_pos) = getInternalWeights(config.internal_weight_initialisation,...
            population(pop_indx).nodes(i),population(pop_indx).nodes(i),config.internal_sparsity);
        
    end
    
    % mask to indicate separation between reservoirs
    population(pop_indx).test_mask{1} = (population(pop_indx).W_scaling == 0);

    % mask to indicate architecture to use. All ones are mutateable
    switch(config.RoR_structure)
        case 'ensemble'
            population(pop_indx).test_mask{2} = false(size(population(pop_indx).test_mask{1}));
        case 'forward_only'
            population(pop_indx).test_mask{2} = triu((population(pop_indx).W_scaling == 0));
        case 'feedback_only'
            population(pop_indx).test_mask{2} = tril((population(pop_indx).W_scaling == 0));
        case 'Graph'
            config.self_loop = [0];               % give node a loop to self. Must be defined as array.
            % node details and connectivity
            t_config = config;
            
            % error checks - currently, all subreservoirs must be equal size
            if contains(config.graph_type,'Lattice')
                if (sqrt(length(config.num_nodes)) ~= round(sqrt(length(config.num_nodes))))
                    error('\n Number of nodes needs to be a square number. \n')
                else
                    t_config.num_nodes = sqrt(size(config.num_nodes,2));
                end
            else
                t_config.num_nodes = size(config.num_nodes,2);
            end

            [t_config,~] = getShape(t_config);
            % find indices for graph weights
            graph_indx = logical(full(adjacency(t_config.G{1})));
                        
            population(pop_indx).test_mask{2} = getArchitectureWeights(graph_indx,config);
        otherwise % or free RoR architecture
            population(pop_indx).test_mask{2} = (population(pop_indx).W_scaling == 0);
    end
    

    % add rand output weights
    if config.add_input_states
        output_units = population(pop_indx).total_units + population(pop_indx).n_input_units;
    else
        output_units = population(pop_indx).total_units;
    end
    
    %apply initialisation technique
    population(pop_indx).output_weights = getInternalWeights(config.output_weight_initialisation...
        ,output_units, population(pop_indx).n_output_units,config.output_connectivity);
    
    % add rand feedback weights
    if config.evolve_feedback_weights
        population(pop_indx).feedback_scaling = 2*rand-1;
        population(pop_indx).feedback_weights = getInternalWeights(config.feedback_weight_initialisation...
            ,population(pop_indx).total_units, population(pop_indx).n_output_units,config.output_connectivity);
    end
    
    % assign blank space for behaviours for CHARC characterisation
    population(pop_indx).behaviours = [];
end
end

function weights = getInternalWeights(weight_initialisation,x_dim,y_dim,sparsity)

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

% currently, all subreservoirs must be equal size
function weights = getArchitectureWeights(graph_indx,config)

weights = zeros(config.total_units);

for n = 1:length(graph_indx)
for i = 1:length(graph_indx)
    for j = 1:length(graph_indx)
        
        if graph_indx(i,j) == 1
            x = (i-1)*config.num_nodes(n)+1:(i)*config.num_nodes(n);
            y = (j-1)*config.num_nodes(n)+1:(j)*config.num_nodes(n);
            
            weights(x,y) = 1;
            
        end
    end
end
end
end