%% This new version of the MM function allows the construction of archiectures of magnets. 
% To build layers of magnets makes sure:
%       layer1 = [res1_nodes, res2_nodes,... ];
%       config.num_nodes = {layer1;layer2; ... etc};
% Or,
%       config.num_nodes = {[25,10,50];[10,5]};
%
% This allows uneven number of reservoirs in each layer, e.g. acting as
% tree network. Alternatively, a pipeline (deep) network can be defined (e.g.
% {[],[],[],..} (mulitple layers) or an ensemble, e.g {[10 10 10 10]} (one layer only)

function population = createMultiMM(config)

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
    
    % add single bias node - this is taken from getAdditionalParameters.m
    population(pop_indx).bias_node = config.bias;
    
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
            
            % add core indx
            population(pop_indx).core_indx = pop_indx;

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
               % geometry details
            if config.evolve_geometry
                                
                [x,y,num_points, area_ratio] = getPolyShape(config.shape_list,config.shape_num_sides,config.rotate_angle,config.ref_point);
                
                % update params
                total_cells_needed = (config.dummy_node_list(layer_indx,res_indx)/area_ratio); 
                square_film_dimensions = ceil(sqrt(total_cells_needed));
                population(pop_indx).layer(layer_indx).nodes(res_indx) = square_film_dimensions.^2;
                config.poly_num = num_points;
                
                population(pop_indx).layer(layer_indx).poly_coord = [x; y]';
                
                
                % % define inputs to sweep through
                [xq,yq] = meshgrid(linspace(0,1,square_film_dimensions),linspace(0,1,square_film_dimensions));
                [in,on] = inpolygon(xq,yq,x,y);
                population(pop_indx).layer(layer_indx).inputs_mask = find(in | on);
                
%                 if config.evolve_poly
%                     population(pop_indx).layer(layer_indx).poly_coord = config.lb + (1-config.lb)*rand(config.poly_num,2);
%                 else
%                     population(pop_indx).layer(layer_indx).geo_width = config.lb + (1-config.lb)*rand;
%                     population(pop_indx).layer(layer_indx).geo_height = config.lb + (1-config.lb)*rand;
%                 end
            else
            population(pop_indx).layer(layer_indx).nodes(res_indx) = config.dummy_node_list(layer_indx,res_indx);%config.total_units_per_layer(res_indx);
        end
            % check length of params
            if length(config.macro_cell_size) > 1
                population(pop_indx).layer(layer_indx).macro_cell_size(res_indx) = config.macro_cell_size(res_indx);
            else
                population(pop_indx).layer(layer_indx).macro_cell_size(res_indx) = config.macro_cell_size;
            end
            if length(config.system_size_z) > 1
                system_size_z = config.system_size_z(res_indx);
            else
                system_size_z = config.system_size_z;
            end
            
            population(pop_indx).layer(layer_indx).system_size(res_indx,1) = (sqrt(population(pop_indx).layer(layer_indx).nodes(res_indx)) * population(pop_indx).layer(layer_indx).macro_cell_size(res_indx))-1;%floor(sqrt(population(pop_indx).nodes(i)* config.macro_cell_size.^3/config.macro_cell_size))-1;
            population(pop_indx).layer(layer_indx).system_size(res_indx,2) = (sqrt(population(pop_indx).layer(layer_indx).nodes(res_indx)) * population(pop_indx).layer(layer_indx).macro_cell_size(res_indx))-1;%floor(sqrt(population(pop_indx).nodes(i)* config.macro_cell_size.^3/config.macro_cell_size))-1;
            population(pop_indx).layer(layer_indx).system_size(res_indx,3) = system_size_z;
            
            % layer thickness
            population(pop_indx).layer(layer_indx).thickness(res_indx) = ceil(rand*10)/10;
            population(pop_indx).layer(layer_indx).minimum_height(res_indx,:) = [0 population(pop_indx).layer(layer_indx).thickness(res_indx)];
            population(pop_indx).layer(layer_indx).maximum_height(res_indx,:) = [population(pop_indx).layer(layer_indx).thickness(res_indx) 1];
                    
        
            % simulation details
            % Input time-period
            population(pop_indx).layer(layer_indx).time_steps_increment = randi([config.time_steps_increment(1) config.time_steps_increment(2)]);%config.time_steps_increment(1);%
            
            %% global params
            population(pop_indx).layer(layer_indx).input_scaling(res_indx)= 2*rand-1;
            population(pop_indx).layer(layer_indx).leak_rate(res_indx) = rand;
            
% assign interpolation time 
            population(pop_indx).layer(layer_indx).interpolation_length = randi([1 config.max_interpolation_length]);
    

            %% Input params
            
            % set positions of magnetic sources. Need maxpos > minpos
            [config.cell_grid_x, config.cell_grid_y] = meshgrid(1:1:(population(pop_indx).layer(layer_indx).system_size(res_indx,1)+1)/population(pop_indx).layer(layer_indx).macro_cell_size(res_indx),1:1:(population(pop_indx).layer(layer_indx).system_size(res_indx,2)+1)/population(pop_indx).layer(layer_indx).macro_cell_size(res_indx));
            population(pop_indx).layer(layer_indx).num_input_loc = population(pop_indx).layer(layer_indx).nodes(res_indx); 
            
            population(pop_indx).layer(layer_indx).xy{res_indx} = [config.cell_grid_x(:) config.cell_grid_y(:)];
            
            % input weights
            if layer_indx == 1
                
                % single reservoir input example
                if config.single_input 
                    input_weights = zeros(population(pop_indx).n_input_units + 1, population(pop_indx).layer(layer_indx).nodes(res_indx));
                    
                    % make rectangle
                    if config.evolve_geometry
                        xy = createRect(population(pop_indx).layer(layer_indx).geo_height,population(pop_indx).layer(layer_indx).geo_width);
                    else
                       xy = createRect(1,1);
                    end
                    % % reset input weights
                    square_film_dimensions = sqrt(population(pop_indx).layer(layer_indx).nodes(res_indx));
                    [xq,yq] = meshgrid(linspace(0,1,square_film_dimensions),linspace(0,1,square_film_dimensions));
                    [in,on] = inpolygon(xq,yq,xy(:,1),xy(:,2));
                    inputs_in_use = find(in | on);

                    input_loc = randperm(length(inputs_in_use),1);
                    weight_value = 2*rand(population(pop_indx).n_input_units + 1,1)-1;
                    
                    input_weights(:,inputs_in_use(input_loc)) = weight_value;
                    population(pop_indx).layer(layer_indx).input_weights{res_indx} = input_weights;
                else
                    population(pop_indx).layer(layer_indx).input_weights{res_indx} =...
                        getWeights(config.input_weight_initialisation,...
                        population(pop_indx).n_input_units + 1,...
                        population(pop_indx).layer(layer_indx).nodes(res_indx),...
                        config.sparsity); 
                end
           
            else
                switch(population(pop_indx).architecture)
                    case 'pipeline'
                        population(pop_indx).layer(layer_indx).input_weights{res_indx} =...
                            getWeights(config.input_weight_initialisation,...
                            config.total_units_per_layer(layer_indx-1)*length(config.read_mag_direction),...
                            population(pop_indx).layer(layer_indx).nodes(res_indx),...
                            config.connecting_sparsity);
                    case {'pipeline_IA','tree'}
                        population(pop_indx).layer(layer_indx).input_weights{res_indx} =...
                            getWeights(config.input_weight_initialisation,...
                            config.total_units_per_layer(layer_indx-1)*length(config.read_mag_direction) + population(pop_indx).n_input_units + config.bias,...
                            population(pop_indx).layer(layer_indx).nodes(res_indx),...
                            config.connecting_sparsity); 
                end
            end
            
                        
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
            if config.multilayer
                population(pop_indx).layer(layer_indx).interfacial_exchange(res_indx) = config.interfacial_exchange(1) + (config.interfacial_exchange(2)-config.interfacial_exchange(1))*rand;
                population(pop_indx).layer(layer_indx).minimum_height(res_indx,:) = [0 population(pop_indx).layer(layer_indx).thickness(res_indx)];
                population(pop_indx).layer(layer_indx).maximum_height(res_indx,:) = [population(pop_indx).layer(layer_indx).thickness(res_indx) 1];
            end
            % boundary params
            population(pop_indx).layer(layer_indx).periodic_boundary(res_indx,:) = config.periodic_boundary;%zeros(1,3);
            %population(pop_indx).layer(layer_indx).periodic_boundary(res_indx,logical(config.periodic_boundary)) = round(rand(1,sum(config.periodic_boundary)));
            
            
            % apply material densities
            if config.evolve_material_density%(res_indx)
                population(pop_indx).layer(layer_indx).material_density(res_indx) = rand.*config.material_density;
            else
                population(pop_indx).layer(layer_indx).material_density(res_indx) = config.material_density;
            end
            
        end
        
    end
    
    %% define connecting matices - will create a large connecting matrix per layer

    % initialise first matrices
%      population(pop_indx).W{1,1} = getWeights(config.internal_weight_initialisation,...
%         population(pop_indx).layer(1).nodes(1),...
%         population(pop_indx).layer(1).nodes(1),...
%         config.internal_sparsity);
    
    % cycle through all other connecting matrices
%     for layer_indx = 1:config.num_layers
%         for res_indx = 1:config.num_res_in_layer
%             for res_indx2 = 1:config.num_res_in_layer
%                 population(pop_indx).W{layer_indx,res_indx} = getWeights(config.internal_weight_initialisation,...
%                     population(pop_indx).layer(layer_indx).nodes(res_indx),...
%                     population(pop_indx).layer(layer_indx).nodes(res_indx2),...
%                     config.internal_sparsity);
%                 
%             end
%         end
%     end
    
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