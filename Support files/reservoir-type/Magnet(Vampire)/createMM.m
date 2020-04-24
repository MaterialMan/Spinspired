function population = createMM(config)

%% Reservoir Parameters
for pop_indx = 1:config.pop_size
    
    % assign details for vampire directories
    population(pop_indx).batch_num = config.test;
    population(pop_indx).core_indx = pop_indx;
    
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
    
    % iterate through subreservoirs
    for i = 1:config.num_reservoirs
        
        % material details
        population(pop_indx).material_type{i} = config.material_type{randi([1 length(config.material_type)])};
        population(pop_indx).material_shape{i} = config.material_shape{randi([1 length(config.material_shape)])};
        switch(population(pop_indx).material_type{i})
            case 'single'
                population(pop_indx).num_materials(i) = 1;
            otherwise
                population(pop_indx).num_materials(i) = 2;
        end
        
        %define num of units
        population(pop_indx).nodes(i) = config.num_nodes(i);
        
        population(pop_indx).system_size(i,1) = (sqrt(population(pop_indx).nodes(i)) * config.macro_cell_size)-1;%floor(sqrt(population(pop_indx).nodes(i)* config.macro_cell_size.^3/config.macro_cell_size))-1;
        population(pop_indx).system_size(i,2) = (sqrt(population(pop_indx).nodes(i)) * config.macro_cell_size)-1;%floor(sqrt(population(pop_indx).nodes(i)* config.macro_cell_size.^3/config.macro_cell_size))-1;
        population(pop_indx).system_size(i,3) = config.system_size_z(i);
        
        % layer thickness
        population(pop_indx).thickness(i) = round(rand*10)/10;
        population(pop_indx).minimum_height(i,:) = [0 population(pop_indx).thickness(i)];
        population(pop_indx).maximum_height(i,:) = [population(pop_indx).thickness(i) 1];
        
        %% global params
        population(pop_indx).input_scaling(i)= 2*rand-1; % not sure about range?
        population(pop_indx).leak_rate(i) = rand;
        
        %% Input params
        % set positions of magnetic sources. Need maxpos > minpos
        [config.cell_grid_x, config.cell_grid_y] = meshgrid(1:1:sqrt(population(pop_indx).nodes(i)));
        population(pop_indx).num_input_loc = population(pop_indx).nodes(i); %population(pop_indx).n_input_units;
        
        population(pop_indx).xy{i} = [config.cell_grid_x(:) config.cell_grid_y(:)];
        
        if config.sparse_input_weights
            input_weights = sprand(population(pop_indx).num_input_loc,  population(pop_indx).n_input_units+1, config.sparsity(i));
            input_weights(input_weights ~= 0) = ...
                2*input_weights(input_weights ~= 0)  - 1;
            
            population(pop_indx).input_weights{i} = input_weights;
        else
            population(pop_indx).input_weights{i}= 2*rand(population(pop_indx).num_input_loc,population(pop_indx).n_input_units+1)-1; % set random field strength
        end
        
        % input widths
        if config.input_widths
            widths = ceil(abs(randn(length(input_weights),1))); %less likely to get big inputs
            widths(widths > round(sqrt(population(pop_indx).nodes(i))/4)) = round(sqrt(population(pop_indx).nodes(i))/4);% cap at 1/6 size of space
        else
            widths = ones(length(input_weights),1);
        end
        population(pop_indx).input_widths{i} = widths; %size of the inputs; pin-point or broad
        
        population(pop_indx).last_state{i} = zeros(1,population(pop_indx).nodes(i));
        
        %% magnet params
        for m = 1: population(pop_indx).num_materials(i)
            population(pop_indx).damping(i,m) = config.damping_parameter(1) + (config.damping_parameter(2)-config.damping_parameter(1))*rand;
            
            population(pop_indx).anisotropy(i,m) = config.anisotropy_parameter(1) + (config.anisotropy_parameter(2)-config.anisotropy_parameter(1))*rand;
            
            population(pop_indx).temperature(i,m) = config.temperature_parameter(1) + (config.temperature_parameter(2)-config.temperature_parameter(1))*rand;
            
            population(pop_indx).exchange(i,m) = config.exchange_parameter(1) + (config.exchange_parameter(2)-config.exchange_parameter(1))*rand;
            
            population(pop_indx).magmoment(i,m) = config.magmoment_parameter(1) + (config.magmoment_parameter(2)-config.magmoment_parameter(1))*rand;    
        end
        population(pop_indx).applied_field_strength(i) = config.applied_field_strength(1) + (config.applied_field_strength(2)-config.applied_field_strength(1))*rand;
        
                
        % random alloy params
        if config.random_alloy(i)
            population(pop_indx).interfacial_exchange(i) = config.exchange_parameter(1) + (config.exchange_parameter(2)-config.exchange_parameter(1))*rand;
            population(pop_indx).alloy_fraction(i) = rand;
        end
        % core shell params
        if config.core_shell(i)
             population(pop_indx).interfacial_exchange(i) = config.exchange_parameter(1) + (config.exchange_parameter(2)-config.exchange_parameter(1))*rand;
             population(pop_indx).shell_size(i,:) = [1 rand];  
             population(pop_indx).particle_size(i) = population(pop_indx).system_size(i,1);
        end

        % apply material densities
        if config.evolve_material_density(i)
            population(pop_indx).material_density(i) = rand.*config.material_density;
        else
            population(pop_indx).material_density = config.material_density;
        end
            
        population(pop_indx).total_units = population(pop_indx).total_units + population(pop_indx).nodes(i);
    end
    
    
    %% define connecting matices
    for i = 1:config.num_reservoirs
        for j= 1:config.num_reservoirs
            
            % initialise connectivity of subrevervoir weights and connecting weights
            if i == j
                population(pop_indx).connectivity(i,j) = config.internal_sparsity;
            else
                population(pop_indx).connectivity(i,j) = config.connecting_sparsity;
            end
            
            switch(config.architecture)
                case 'ensemble'
                    
                    % no connectivity
                    population(pop_indx).connectivity(i,j) = 0;
                    population(pop_indx).W_scaling(i,j) = 0;
                    population(pop_indx).W{i,j} = zeros(population(pop_indx).nodes(i), population(pop_indx).nodes(j));
                    
                case {'pipeline','pipeline_IA'}  % pipeline structure
                    if  j == i+1
                        switch(config.internal_weight_initialisation)
                            case 'norm' % normal distribution
                                internal_weights = sprandn(population(pop_indx).nodes(i), population(pop_indx).nodes(j), population(pop_indx).connectivity(i,j));
                                
                            case 'uniform' % uniform dist between -1 and 1
                                internal_weights = sprand(population(pop_indx).nodes(i), population(pop_indx).nodes(j), population(pop_indx).connectivity(i,j));
                                internal_weights(internal_weights ~= 0) = ...
                                    2*internal_weights(internal_weights ~= 0)  - 1;
                        end
                        
                        % assign scaling for inner weights
                        population(pop_indx).W_scaling(i,j) = 2*rand-1;
                        population(pop_indx).W{i,j} = internal_weights;
                    else
                        % if self, do nothing
                        population(pop_indx).connectivity(i,j) = 0;
                        population(pop_indx).W_scaling(i,j) = 0;
                        population(pop_indx).W{i,j} = zeros(population(pop_indx).nodes(i), population(pop_indx).nodes(j));
                    end
                    
            end
        end
        
    end
    
    % add rand output weights
    if config.add_input_states
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units + population(pop_indx).n_input_units, population(pop_indx).n_output_units)-1;
    else
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units, population(pop_indx).n_output_units)-1;
    end
    
    population(pop_indx).behaviours = [];
    
end
