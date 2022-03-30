function population = createSTO(config)
% To build layers of magnets makes sure:
%       layer1 = [res1_nodes, res2_nodes,... ];
%       config.num_nodes = {layer1;layer2; ... etc};
% Or,
%       config.num_nodes = {[25,10,50];[10,5]};
%
% This allows uneven number of reservoirs in each layer, e.g. acting as
% tree network. Alternatively, a pipeline (deep) network can be defined (e.g.
% {[],[],[],..} (mulitple layers) or an ensemble, e.g {[10 10 10 10]} (one layer only)


%set up the new batch path
config.batch_path = getMMPath(config.test);

%% Reservoir Parameters
for pop_indx = 1:config.pop_size
    
    % assign details for vampire directories
    population(pop_indx).batch_num = config.test;
    population(pop_indx).batch_path = config.batch_path;
    
    % add performance records
    population(pop_indx).train_error = 1;
    population(pop_indx).val_error = 1;
    population(pop_indx).test_error = 1;
    

    % assign input/output count
    if isempty(config.train_input_sequence)
        population(pop_indx).n_input_units = 1;
        population(pop_indx).n_output_units = 1;
    else
        population(pop_indx).n_input_units = config.task_num_inputs;
        population(pop_indx).n_output_units = config.task_num_outputs;
    end
    
    %track total nodes in use
    population(pop_indx).total_units = config.total_units;
    
    % assign architecture and material details
    population(pop_indx).architecture = config.architecture;
    
    % iterate through layers
    for layer_indx = 1:config.num_layers
        
        % iterate through subreservoirs in layer
        for res_indx = 1:config.num_res_in_layer(layer_indx)
            
            population(pop_indx).core_indx = layer_indx;
            
            %define num of units
            population(pop_indx).layer(layer_indx).nodes(res_indx) = config.dummy_node_list(layer_indx,res_indx);%config.total_units_per_layer(res_indx);
            
            % add single bias node - this is taken from getAdditionalParameters.m
            population(pop_indx).layer(layer_indx).bias_node(res_indx) = config.bias;
            
            population(pop_indx).layer(layer_indx).unit_cell_size(res_indx,:) = repmat((population(pop_indx).layer(layer_indx).nodes(res_indx)*16),1,3);
            population(pop_indx).layer(layer_indx).system_size(res_indx,:) = (population(pop_indx).layer(layer_indx).unit_cell_size(res_indx,:)/10)-1;
            
            % simulation details
            % Input time-period
            population(pop_indx).layer(layer_indx).time_steps_increment(res_indx) = randi([config.time_steps_increment(1) config.time_steps_increment(2)]);
            population(pop_indx).layer(layer_indx).interpolation_length(res_indx) = randi([config.interpolation_length(1) config.interpolation_length(2)]);
            
            %% global params
            population(pop_indx).layer(layer_indx).input_scaling(res_indx)= 2*rand-1;
            population(pop_indx).layer(layer_indx).leak_rate(res_indx) = rand;
            
            %% Input params
            % define array layout
            population(pop_indx).layer(layer_indx).array_mask{res_indx} = full(adjacency(config.G{layer_indx}{res_indx}));
            
            % record connectivity of array
            [x,y] = find(adjacency(config.G{layer_indx}{res_indx}));
            population(pop_indx).layer(layer_indx).xy{res_indx} = [x,y]-1;
            
            % input weights
            if layer_indx == 1
                % inputweights
                population(pop_indx).layer(layer_indx).input_weights{res_indx} =...
                    getWeights(config.input_weight_initialisation,...
                    ones(population(pop_indx).n_input_units + 1,population(pop_indx).layer(layer_indx).nodes(res_indx)),...
                    population(pop_indx).n_input_units + 1,...
                    population(pop_indx).layer(layer_indx).nodes(res_indx),...
                    config.sparsity);
            else
                switch(population(pop_indx).architecture)
                    case 'pipeline'
                        population(pop_indx).layer(layer_indx).input_weights{res_indx} =...
                            getWeights(config.input_weight_initialisation,...
                            ones(config.total_units_per_layer(layer_indx-1)*length(config.read_mag_direction),population(pop_indx).layer(layer_indx).nodes(res_indx)),...
                            config.total_units_per_layer(layer_indx-1)*length(config.read_mag_direction),...
                            population(pop_indx).layer(layer_indx).nodes(res_indx),...
                            config.connecting_sparsity);
                    case {'pipeline_IA','tree'}
                        population(pop_indx).layer(layer_indx).input_weights{res_indx} =...
                            getWeights(config.input_weight_initialisation,...
                            ones(config.total_units_per_layer(layer_indx-1)*length(config.read_mag_direction)+ population(pop_indx).n_input_units + config.bias,population(pop_indx).layer(layer_indx).nodes(res_indx)),...
                            config.total_units_per_layer(layer_indx-1)*length(config.read_mag_direction) + population(pop_indx).n_input_units + config.bias,...
                            population(pop_indx).layer(layer_indx).nodes(res_indx),...
                            config.connecting_sparsity);
                end
            end
            
            population(pop_indx).layer(layer_indx).last_state{res_indx} = zeros(1,population(pop_indx).layer(layer_indx).nodes(res_indx));
            
            population(pop_indx).layer(layer_indx).num_materials(res_indx) = 1;
            
            %% magnet params
            for m = 1: population(pop_indx).layer(layer_indx).num_materials(res_indx)
                population(pop_indx).layer(layer_indx).damping(res_indx,m) = config.damping_parameter(1) + (config.damping_parameter(2)-config.damping_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).anisotropy(res_indx,m) = config.anisotropy_parameter(1) + (config.anisotropy_parameter(2)-config.anisotropy_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).temperature(res_indx,m) = config.temperature_parameter(1) + (config.temperature_parameter(2)-config.temperature_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).exchange(res_indx,m) = config.exchange_parameter(1) + (config.exchange_parameter(2)-config.exchange_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).magmoment(res_indx,m) = config.magmoment_parameter(1) + (config.magmoment_parameter(2)-config.magmoment_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).applied_field_strength(res_indx,m) = config.applied_field_strength(1) + (config.applied_field_strength(2)-config.applied_field_strength(1))*rand;
                
                population(pop_indx).layer(layer_indx).spin_transfer_relaxation_torque(res_indx,m) = config.spin_transfer_relaxation_torque(1) + (config.spin_transfer_relaxation_torque(2)-config.spin_transfer_relaxation_torque(1))*rand;
                
                population(pop_indx).layer(layer_indx).spin_transfer_precession_torque(res_indx,m) = config.spin_transfer_precession_torque(1) + (config.spin_transfer_precession_torque(2)-config.spin_transfer_precession_torque(1))*rand;
                
                population(pop_indx).layer(layer_indx).spin_transfer_torque_asymmetry(res_indx,m) = config.spin_transfer_torque_asymmetry(1) + (config.spin_transfer_torque_asymmetry(2)-config.spin_transfer_torque_asymmetry(1))*rand;
                
            end
            
            population(pop_indx).layer(layer_indx).W_scaling(res_indx) = 2*rand-1;
            
            if config.hetero_exchange_interaction
                weight_initialisation = 'norm';
            else
                weight_initialisation = 'ones';
            end
            population(pop_indx).layer(layer_indx).W{res_indx} =...
                getWeights(weight_initialisation,...
                population(pop_indx).layer(layer_indx).array_mask{res_indx},...
                population(pop_indx).layer(layer_indx).nodes(res_indx),...
                population(pop_indx).layer(layer_indx).nodes(res_indx),...
                1);
            
            % remove node 1 connections
            population(pop_indx).layer(layer_indx).W{res_indx}(1,:) = 0;
            population(pop_indx).layer(layer_indx).W{res_indx}(:,1) = 0;
        end
    end
    
    
    % add rand output weights
    if config.add_input_states
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units*length(config.read_mag_direction) + population(pop_indx).n_input_units, population(pop_indx).n_output_units)-1;
    else
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units*length(config.read_mag_direction), population(pop_indx).n_output_units)-1;
    end
    
    population(pop_indx).behaviours = [];
    
end

end


function weights = getWeights(weight_initialisation,array_mask,x_dim,y_dim,sparsity)

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
    case 'full'
        weights = rand(x_dim, y_dim);
    case 'ones'
        weights = ones(x_dim, y_dim);
end

% rescale - between 0 and 1 ('rescale', [0 1]) for inner weights and -1 and for
% connecting weights ('rescale', [-1 1])
% tmp_config.preprocess = scale;
% tmp_config.preprocess_shift = shift;
% weights = featureNormailse(weights,tmp_config);

% if not defined in connectivity mask, delete
weights(~array_mask) = 0;

end

