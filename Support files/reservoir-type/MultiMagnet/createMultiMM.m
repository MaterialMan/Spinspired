%% This new version of the MM function allows the construction of archiectures of magnets. 
% To build layers of magnetics makes sure:
%       layer1 = [res1_nodes, res2_nodes,... ];
%       config.num_nodes = {layer1;layer2; ... etc};
% Or,
%       config.num_nodes = {[25,10,50];[10,5]};

function population = createMultiMM(config)

%% Reservoir Parameters
for pop_indx = 1:config.pop_size
    
    % assign details for vampire directories
    population(pop_indx).batch_num = config.test;
    
    % add performance records
    population(pop_indx).train_error = 1;
    population(pop_indx).val_error = 1;
    population(pop_indx).test_error = 1;
    
    % add single bias node - this is taken from getAdditionalParameters.m
    population(pop_indx).bias_node = config.bias;
    
    % assign input/output count
    if isempty(config.train_input_sequence)
        population(pop_indx).n_input_units = 1;
        population(pop_indx).n_output_units = 1;
    else
        population(pop_indx).n_input_units = size(config.train_input_sequence,2);
        population(pop_indx).n_output_units = size(config.train_output_sequence,2);
    end
    
    %track total nodes in use
    population(pop_indx).total_units = 0;
    
    % assign architecture and material details
    population(pop_indx).architecture = config.architecture;
    
    % iterate through layers
    for layer_indx = 1:config.num_layers
        
        % iterate through subreservoirs
        for res_indx = 1:config.num_reservoirs(layer_indx)
            
            if config.num_reservoirs(layer_indx) > 1
                population(pop_indx).layer(layer_indx).core_indx(res_indx) = res_indx;
            else
                population(pop_indx).layer(layer_indx).core_indx = pop_indx;
            end
            
            % material details
            population(pop_indx).layer(layer_indx).material_type{res_indx} = config.material_type{randi([1 length(config.material_type)])};
            population(pop_indx).layer(layer_indx).material_shape{res_indx} = config.material_shape{randi([1 length(config.material_shape)])};
            switch(population(pop_indx).layer(layer_indx).material_type{res_indx})
                case {'multilayer','core_shell', 'random_alloy'}
                    population(pop_indx).layer(layer_indx).num_materials(res_indx) = 2;
                otherwise
                    population(pop_indx).layer(layer_indx).num_materials(res_indx) = 1;
            end
            
            %define num of units
            population(pop_indx).layer(layer_indx).nodes(res_indx) = config.num_nodes{layer_indx}(res_indx);
            
            % check length of params
            if length(config.macro_cell_size) > 1
                population(pop_indx).layer(layer_indx).macro_cell_size(res_indx) = config.macro_cell_size(res_indx);
            else
                population(pop_indx).layer(layer_indx).macro_cell_size(res_indx) = config.macro_cell_size;
            end
            if (config.system_size_z) > 1
                system_size_z = config.system_size_z(res_indx);
            else
                system_size_z = config.system_size_z;
            end
            
            population(pop_indx).layer(layer_indx).system_size(res_indx,1) = (sqrt(population(pop_indx).layer(layer_indx).nodes(res_indx)) * population(pop_indx).layer(layer_indx).macro_cell_size(res_indx))-1;%floor(sqrt(population(pop_indx).nodes(i)* config.macro_cell_size.^3/config.macro_cell_size))-1;
            population(pop_indx).layer(layer_indx).system_size(res_indx,2) = (sqrt(population(pop_indx).layer(layer_indx).nodes(res_indx)) * population(pop_indx).layer(layer_indx).macro_cell_size(res_indx))-1;%floor(sqrt(population(pop_indx).nodes(i)* config.macro_cell_size.^3/config.macro_cell_size))-1;
            population(pop_indx).layer(layer_indx).system_size(res_indx,3) = system_size_z;
            
            % layer thickness
            population(pop_indx).layer(layer_indx).thickness(res_indx) = round(rand*10)/10;
            population(pop_indx).layer(layer_indx).minimum_height(res_indx,:) = [0 population(pop_indx).layer(layer_indx).thickness(res_indx)];
            population(pop_indx).layer(layer_indx).maximum_height(res_indx,:) = [population(pop_indx).layer(layer_indx).thickness(res_indx) 1];
            
            % simulation details
            % Input time-period
            population(pop_indx).layer(layer_indx).time_steps_increment = randi([config.time_steps_increment(1) config.time_steps_increment(2)]);%config.time_steps_increment(1);%
            
            %% global params
            population(pop_indx).layer(layer_indx).input_scaling(res_indx)= 2*config.input_scaler*rand - config.input_scaler;%*(2*rand-1); % not sure about range?
            population(pop_indx).layer(layer_indx).leak_rate(res_indx) = rand;
            
            %% Input params
            
            % set positions of magnetic sources. Need maxpos > minpos
            [config.cell_grid_x, config.cell_grid_y] = meshgrid(1:1:(population(pop_indx).layer(layer_indx).system_size(res_indx,1)+1)/population(pop_indx).layer(layer_indx).macro_cell_size(res_indx),1:1:(population(pop_indx).layer(layer_indx).system_size(res_indx,2)+1)/population(pop_indx).layer(layer_indx).macro_cell_size(res_indx));
            population(pop_indx).layer(layer_indx).num_input_loc = population(pop_indx).layer(layer_indx).nodes(res_indx); %population(pop_indx).n_input_units;
            
            population(pop_indx).layer(layer_indx).xy{res_indx} = [config.cell_grid_x(:) config.cell_grid_y(:)];
            
            
            % inputweights
            population(pop_indx).layer(layer_indx).input_weights{res_indx} =...
                getWeights(config.input_weight_initialisation,...
                population(pop_indx).layer(layer_indx).nodes(res_indx),...
                population(pop_indx).n_input_units + config.bias,...
                config.sparsity);
                        
            % input widths
            if config.input_widths
                widths = ceil(abs(randn(length(population(pop_indx).layer(layer_indx).input_weights{res_indx}),1))); %less likely to get big inputs
                widths(widths > round(sqrt(population(pop_indx).layer(layer_indx).nodes(res_indx))/4)) = round(sqrt(population(pop_indx).layer(layer_indx).nodes(res_indx))/4);% cap at 1/6 size of space
            else
                widths = ones(length(population(pop_indx).layer(layer_indx).input_weights{res_indx}),1);
            end
            population(pop_indx).layer(layer_indx).input_widths{res_indx} = widths; %size of the inputs; pin-point or broad
            
            population(pop_indx).layer(layer_indx).last_state{res_indx} = zeros(1,population(pop_indx).layer(layer_indx).nodes(res_indx));
            
            %% magnet params
            for m = 1: population(pop_indx).layer(layer_indx).num_materials(res_indx)
                population(pop_indx).layer(layer_indx).damping(res_indx,m) = config.damping_parameter(1) + (config.damping_parameter(2)-config.damping_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).anisotropy(res_indx,m) = config.anisotropy_parameter(1) + (config.anisotropy_parameter(2)-config.anisotropy_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).temperature(res_indx,m) = config.temperature_parameter(1) + (config.temperature_parameter(2)-config.temperature_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).temperature_rescaling_exponent(res_indx,m) = config.temperature_rescaling_exponent(1) + (config.temperature_rescaling_exponent(2)-config.temperature_rescaling_exponent(1))*rand;
                
                population(pop_indx).layer(layer_indx).temperature_rescaling_curie_temperature(res_indx,m) = config.temperature_rescaling_curie_temperature(1) + (config.temperature_rescaling_curie_temperature(2)-config.temperature_rescaling_curie_temperature(1))*rand;
                
                population(pop_indx).layer(layer_indx).exchange(res_indx,m) = config.exchange_parameter(1) + (config.exchange_parameter(2)-config.exchange_parameter(1))*rand;
                
                population(pop_indx).layer(layer_indx).magmoment(res_indx,m) = config.magmoment_parameter(1) + (config.magmoment_parameter(2)-config.magmoment_parameter(1))*rand;
            end
            population(pop_indx).layer(layer_indx).applied_field_strength(res_indx) = config.applied_field_strength(1) + (config.applied_field_strength(2)-config.applied_field_strength(1))*rand;
            
            
            % random alloy params
            if config.random_alloy%(res_indx)
                population(pop_indx).layer(layer_indx).interfacial_exchange(res_indx) = config.interfacial_exchange(1) + (config.interfacial_exchange(2)-config.interfacial_exchange(1))*rand;
                population(pop_indx).layer(layer_indx).alloy_fraction(res_indx) = rand;
            end
            % core shell params
            if config.core_shell%(res_indx)
                population(pop_indx).layer(layer_indx).interfacial_exchange(res_indx) = config.interfacial_exchange(1) + (config.interfacial_exchange(2)-config.interfacial_exchange(1))*rand;
                population(pop_indx).layer(layer_indx).shell_size(res_indx,:) = [1 rand];
                population(pop_indx).layer(layer_indx).particle_size(res_indx) = population(pop_indx).layer(layer_indx).system_size(res_indx,1);
            end
            
            % boundary params
            population(pop_indx).layer(layer_indx).periodic_boundary(res_indx,:) = zeros(1,3);
            population(pop_indx).layer(layer_indx).periodic_boundary(res_indx,logical(config.periodic_boundary)) = round(rand(1,sum(config.periodic_boundary)));
            
            
            % apply material densities
            if config.evolve_material_density%(res_indx)
                population(pop_indx).layer(layer_indx).material_density(res_indx) = rand.*config.material_density;
            else
                population(pop_indx).layer(layer_indx).material_density(res_indx) = config.material_density;
            end
            
            population(pop_indx).total_units = population(pop_indx).total_units + population(pop_indx).layer(layer_indx).nodes(res_indx);
        end
        
    end
    
    %% define connecting matices - will create a large connecting matrix per layer

    % initialise first matrices
%      population(pop_indx).W{1,1} = getWeights(config.internal_weight_initialisation,...
%         population(pop_indx).layer(1).nodes(1),...
%         population(pop_indx).layer(1).nodes(1),...
%         config.internal_sparsity);
    
    % cycle through all other connecting matrices
    for layer_indx = 1:config.num_layers
        for res_indx = 1:length(population(pop_indx).layer(layer_indx).nodes)
            for res_indx2 = 1:length(population(pop_indx).layer(layer_indx).nodes)
                population(pop_indx).W{layer_indx,res_indx} = getWeights(config.internal_weight_initialisation,...
                    population(pop_indx).layer(layer_indx).nodes(res_indx),...
                    population(pop_indx).layer(layer_indx).nodes(res_indx2),...
                    config.internal_sparsity);
                
            end
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

function weights = getWeights(weight_initialisation,x_dim,y_dim,sparsity)

switch(weight_initialisation)
    case 'norm' % normal distribution
        weights = sprandn(x_dim,y_dim,sparsity);
    case 'uniform' % uniform dist between -1 and 1
       weights = sprand(x_dim,y_dim,sparsity);   
       weights(weights ~= 0) = ...
        2*weights(weights ~= 0)  - 1;
    case 'orth'
        weights = orth(full(sprand(x_dim,y_dim,sparsity)));
end
end