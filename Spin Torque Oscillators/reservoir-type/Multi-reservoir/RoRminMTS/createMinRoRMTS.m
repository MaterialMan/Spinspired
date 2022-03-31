function population = createMinRoRMTS(config)


%% Reservoir Parameters
for pop_indx = 1:config.pop_size
    
    % add performance records
    population(pop_indx).train_error = 1;
    population(pop_indx).val_error = 1;
    population(pop_indx).test_error = 1;
    
    % add single bias node
    population(pop_indx).bias_node = config.bias_node;
    
    % assign input/output count
    if isempty(config.train_input_sequence)
        population(pop_indx).n_input_units = 1;
        population(pop_indx).n_output_units = 1;
    else
        population(pop_indx).n_input_units = config.task_num_inputs;
        population(pop_indx).n_output_units = config.task_num_outputs;
    end
    
    
    %define num of units
    population(pop_indx).nodes = config.num_nodes;
    population(pop_indx).total_units = sum(config.num_nodes);
    
    % Scaling and leak rate
    population(pop_indx).init_input_scaling = 2*rand(1,length(config.num_nodes))-1; % can expand to have for each input and bias: population(pop_indx).n_input_units+1
    population(pop_indx).init_W_scaling = 2*rand(1,length(config.num_nodes))-1;
    population(pop_indx).init_leak_rate = rand(1,sum(config.num_nodes));
    
    % update cycle per node
    population(pop_indx).init_update_cycle = randi([0 config.max_update_cycle],1,size(config.num_nodes,2));
    population(pop_indx).init_max_update_cycle = config.max_update_cycle;
    
    % input delay
    population(pop_indx).init_input_delay = randi([0 config.max_input_delay],1,size(config.num_nodes,2));
 
    % assign interpolation time 
    population(pop_indx).interpolation_length = randi([1 config.max_interpolation_length]); %,1,size(config.num_nodes,2)
    
    % quantization
    if config.quantized_state(1) > 0
        population(pop_indx).m = randi([config.quantized_state(1) config.quantized_state(2)],1,size(config.num_nodes,2));
    end
    
    %inputweights
    population(pop_indx).input_weights = getInternalWeights(population(pop_indx),config.input_weight_initialisation...
        ,population(pop_indx).n_input_units+1, population(pop_indx).total_units,config.sparsity,config);
    
    %assign different activations, if necessary
    if config.multi_activ
        population(pop_indx).activ_Fcn_indx = randi(length(config.activ_list),1,population(pop_indx).total_units);
    else
        population(pop_indx).activ_Fcn_indx = randi(length(config.activ_list));
    end
    
    population(pop_indx).last_state = zeros(50,population(pop_indx).total_units);
    
    %% weights and connectivity of all reservoirs
    % end of first subres
    end_pos = population(pop_indx).nodes(1);
    
    % intialise W and W_scale
    population(pop_indx).W_scaling(1:end_pos,1:end_pos) = population(pop_indx).init_W_scaling(1);
    population(pop_indx).input_scaling(:,1:end_pos) = population(pop_indx).init_input_scaling(:,1);
    
    if config.multi_leak_rate 
        population(pop_indx).leak_rate(1:end_pos) = rand(1,sum(config.num_nodes)); % assign indvidual leak rates per node
    else
        population(pop_indx).leak_rate(1:end_pos) = population(pop_indx).init_leak_rate(1);
    end
    
    % initialise update cycles
    if config.per_node_time_scale
        population(pop_indx).init_update_cycle = ones(1,sum(config.num_nodes));%randi([1 config.max_update_cycle],1,sum(config.num_nodes));% assign indvidual update rates per node
        population(pop_indx).update_cycle = population(pop_indx).init_update_cycle;
        %population(pop_indx).max_update_cycle = population(pop_indx).init_max_update_cycle;
    else
        population(pop_indx).update_cycle(1:end_pos) = population(pop_indx).init_update_cycle(1); % need to test
        %population(pop_indx).max_update_cycle(1:end_pos) = population(pop_indx).init_max_update_cycle(1);
    end
    population(pop_indx).max_update_cycle = population(pop_indx).init_max_update_cycle;
    
    % input delay
    if config.per_input_delay
        population(pop_indx).init_input_delay = ones(1,sum(config.num_nodes));
        population(pop_indx).input_delay = population(pop_indx).init_input_delay;
    else
        population(pop_indx).input_delay(1:end_pos) = population(pop_indx).init_input_delay(1); % need to test
    end
    population(pop_indx).max_input_delay = config.max_input_delay;
    
    if config.dipole_fields        
        population(pop_indx).perim_value = [-config.W_scaling, -config.self_loop_scaling, config.corner_scale];
        population(pop_indx).self_feed_back_loop = config.self_feed_back_loop;
        population(pop_indx).dipole_scale = config.dipole_scale;
        population(pop_indx).num_dipoles = config.num_dipoles;
    end
    
    % assign first internal weights
    population(pop_indx).W = getInternalWeights(population(pop_indx),config.internal_weight_initialisation,...
        population(pop_indx).nodes(1),population(pop_indx).nodes(1),config.internal_sparsity,config,1);
       
    if contains(config.internal_weight_initialisation,["graph","weight_fcn","dipole_fields"])
        population(pop_indx).W_scaling(1:end_pos,1:end_pos) = adjacency(config.G{1}).*population(pop_indx).init_W_scaling(1);
    else
        population(pop_indx).W_scaling(1:end_pos,1:end_pos) = population(pop_indx).init_W_scaling(1);
    end

    % cycle through subreservoirs
    for i = 2:length(population(pop_indx).nodes)
        
        start_pos = end_pos+1;
        end_pos = start_pos + population(pop_indx).nodes(i)-1;
        
        population(pop_indx).input_scaling(:,start_pos:end_pos) = population(pop_indx).init_input_scaling(:,i);
        %population(pop_indx).W_scaling(start_pos:end_pos,start_pos:end_pos) = population(pop_indx).init_W_scaling(i);
        
        if ~config.multi_leak_rate
            population(pop_indx).leak_rate(start_pos:end_pos) = population(pop_indx).init_leak_rate(i);
        end
        
        if ~config.per_node_time_scale
            population(pop_indx).update_cycle(start_pos:end_pos) = population(pop_indx).init_update_cycle(i);
%            population(pop_indx).max_update_cycle(start_pos:end_pos) = population(pop_indx).init_max_update_cycle(i);
        end
        
        if ~config.per_input_delay
            population(pop_indx).input_delay(start_pos:end_pos) = population(pop_indx).init_input_delay(i);
        end
        
        % assign new subres weights
        population(pop_indx).W(start_pos:end_pos,start_pos:end_pos) = getInternalWeights(population(pop_indx),config.internal_weight_initialisation,...
            population(pop_indx).nodes(i),population(pop_indx).nodes(i),config.internal_sparsity,config,i);
        
        % assign scaling
        if contains(config.internal_weight_initialisation,["graph","weight_fcn","dipole_fields"])
            population(pop_indx).W_scaling(start_pos:end_pos,start_pos:end_pos) = adjacency(config.G{1}).*population(pop_indx).init_W_scaling(i);
        else
            population(pop_indx).W_scaling(start_pos:end_pos,start_pos:end_pos) = population(pop_indx).init_W_scaling(i);
        end

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
            %config.self_loop = [0];               % give node a loop to self. Must be defined as array.
            % node details and connectivity
            t_config = config;
            t_config.graph_type = config.structure_type;
            t_config.self_loop = [0]; 
            
            % error checks - currently, all subreservoirs must be equal size
            if contains(t_config.graph_type,'Lattice')
                if (sqrt(length(t_config.num_nodes)) ~= round(sqrt(length(t_config.num_nodes))))
                    error('\n Number of nodes needs to be a square number. \n')
                else
                    t_config.num_nodes = sqrt(size(t_config.num_nodes,2));
                end
            else
                t_config.num_nodes = size(t_config.num_nodes,2);
            end

            [t_config,~] = getShape(t_config);
            % find indices for graph weights
            graph_indx = logical(full(adjacency(t_config.G{1})));
                        
            population(pop_indx).test_mask{2} = getArchitectureWeights(graph_indx,config);
            
        otherwise % or free RoR architecture
             if contains(config.internal_weight_initialisation,["graph","weight_fcn","dipole_fields"])
                population(pop_indx).test_mask{2} = (population(pop_indx).W_scaling ~= 0);
            else
                population(pop_indx).test_mask{2} = (population(pop_indx).W_scaling == 0);
            end

    end
    

    % add rand output weights
    if config.add_input_states
        output_units = population(pop_indx).total_units + population(pop_indx).n_input_units;
    else
        output_units = population(pop_indx).total_units;
    end
    
    %apply initialisation technique
    population(pop_indx).output_weights = getInternalWeights(population(pop_indx),config.output_weight_initialisation...
        ,output_units, population(pop_indx).n_output_units,config.output_connectivity,config);
    
    % add rand feedback weights
    if config.evolve_feedback_weights
        population(pop_indx).feedback_scaling = 2*rand-1;
        population(pop_indx).feedback_weights = getInternalWeights(population(pop_indx),config.feedback_weight_initialisation...
            ,population(pop_indx).total_units, population(pop_indx).n_output_units,config.output_connectivity,config);
    end
    
    % assign blank space for behaviours for CHARC characterisation
    population(pop_indx).behaviours = [];
end
end

function weights = getInternalWeights(individual,weight_initialisation,x_dim,y_dim,sparsity,config,indx)

if nargin < 7
    indx = 1;
end

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
        
        case 'graph'
        
        weights = adjacency(config.G{indx}).*(2.*rand(size(adjacency(config.G{indx})))-1);

    case 'weight_fcn'
        
        graph_indx = logical(full(adjacency(config.G{1})));
        
        [x,y] = meshgrid(1:config.num_nodes);
        
        f = ones(config.num_nodes);%config.weight_fcn(x,y);
        
        f = f./max(max(f));
        
        % set feedback weights
        idx = logical(eye(config.num_nodes,config.num_nodes));
        f(idx) = config.self_loop;

        % out nodes
        sq = sqrt(config.num_nodes);
        perim_nodes = [1:sq; 1:sq:config.num_nodes; sq:sq:config.num_nodes;
            (1:sq)+(config.num_nodes-sq)];
        %f(perim_nodes,perim_nodes) = config.perim_value;
        
        tmp_fperim = f(perim_nodes,perim_nodes);
        % define perimeter self-connections
        tmp_fperim(logical(eye(size(tmp_fperim,1)))) = diag(tmp_fperim)*config.perim_value(1);
        % set connections to nearest-neighbours
        A = ones(size(tmp_fperim,1));
        B1 = triu(A,1);
        B2 = tril(A,-1);
        tmp_fperim(logical(B1)) = config.perim_value(2);
        tmp_fperim(logical(B2)) = config.perim_value(3);
        
        % put back into perimeter nodes
        f(perim_nodes,perim_nodes) = tmp_fperim;
        
        weights = f.*graph_indx;
        
        % plot weights
        subplot(1,2,1)
        imagesc(f)
        subplot(1,2,2)
        g = digraph(weights);
        p = plot(g);
        g.Edges.EdgeColors = g.Edges.Weight;
        p.EdgeCData = g.Edges.EdgeColors;
        drawnow
        
    case 'dipole_fields'
                
        weights = adjacency(config.G{indx});
        
        % set strength of dipole fields
        dipole_scale = individual.dipole_scale; %0.2
        dipole_field_value = (1./((1:config.num_nodes(indx)).^3))*dipole_scale;%1./exp((1:config.config.num_nodes(indx))*1.1)*3;% 1./exp((1:config.config.num_nodes(indx))*1.1)*3
%         for d = 2:2:individual.num_dipoles
%              dipole_field_value(d) = -dipole_field_value(d);
%          end
               
        for n = 1:config.num_nodes(indx)
            source_node = n;
            for i = 1:individual.num_dipoles
                [nn,dist] = nearest(config.G{indx},source_node,i,'Method','unweighted');
                weights(sub2ind(size(weights),repmat(source_node,length(nn(dist == i)),1),nn(dist == i))) = dipole_field_value(i);
                weights(sub2ind(size(weights),nn(dist == i),repmat(source_node,length(nn(dist == i)),1))) = dipole_field_value(i);
            end
        end
                
        % set feedback weights
%         idx = logical(eye(config.num_nodes(indx),config.num_nodes(indx)));
%         weights(idx) = individual.self_feed_back_loop;

%         % update corner links - rescale corner connections
%         for i = 1:config.num_nodes(indx)
%             
%             if i-sqrt(config.num_nodes(indx))-1 > 0
%                 individual.W_scaling(i,i-sqrt(config.num_nodes(indx))-1) = individual.W_scaling(i,i-sqrt(config.num_nodes(indx))-1)*config.corner_scale;
%                 
%             end
%             if i-sqrt(config.num_nodes(indx))+1 > 0
%                 individual.W_scaling(i,i-sqrt(config.num_nodes(indx))+1) = individual.W_scaling(i,i-sqrt(config.num_nodes(indx))+1)*config.corner_scale;
%                 
%             end
%             if i+sqrt(config.num_nodes(indx))-1 > 0 && i > 1 && i+sqrt(config.num_nodes(indx))+1 < config.num_nodes(indx)
%                 individual.W_scaling(i,i+sqrt(config.num_nodes(indx))-1) = individual.W_scaling(i,i+sqrt(config.num_nodes(indx))-1)*config.corner_scale;
%                 
%             end
%             if i+sqrt(config.num_nodes(indx))+1 > 0 && i+sqrt(config.num_nodes(indx))+1 < config.num_nodes(indx)
%                 individual.W_scaling(i,i+sqrt(config.num_nodes(indx))+1) = individual.W_scaling(i,i+sqrt(config.num_nodes(indx))+1)*config.corner_scale;
%                 
%             end
%         end
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