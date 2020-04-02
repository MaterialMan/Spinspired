%% create_ReservoirName_.m
% Template function to define reservoir parameters. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
% Tip: Try maintain ordering as structs keep this ordering.

% This is called by the @config.createFcn pointer.

function population = createOregonator(config)

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
        population(pop_indx).n_input_units = size(config.train_input_sequence,2);
        population(pop_indx).n_output_units = size(config.train_output_sequence,2);
    end
    
    % iterate through subreservoirs
    for i = 1:config.num_reservoirs
        
        %define num of units
        population(pop_indx).nodes(i) = config.num_nodes(i);
        
        % oregonator params - add anything that needs to be evolved
        population(pop_indx).vesicle_grid_size(i) = sqrt(config.num_nodes(i)); 
        population(pop_indx).height(i) = 2*config.vesicle_radius*population(pop_indx).vesicle_grid_size(i);
        population(pop_indx).width(i) = 2*config.vesicle_radius*population(pop_indx).vesicle_grid_size(i);               
        
        population(pop_indx).time_period(i) = randi([config.max_time_period config.max_time_period]);
        population(pop_indx).input_length(i) = randi([1 population(pop_indx).time_period(i)]);
        population(pop_indx).input_widths(i) = randi([1 config.vesicle_radius-1]);
        
        % Scaling and leak rate
        population(pop_indx).input_scaling(i) = 0.05*rand; %increases nonlinearity
        population(pop_indx).leak_rate(i) = rand;
        
        % input weights
        switch(config.input_weight_initialisation)
            case 'norm' % normal distribution
                if config.sparse_input_weights
                    input_weights = sprandn(population(pop_indx).nodes(i),  population(pop_indx).n_input_units+1, config.sparsity(i));
                else
                    input_weights = randn(population(pop_indx).nodes(i),  population(pop_indx).n_input_units+1);
                end
            case 'uniform' % uniform dist between -1 and 1
                if config.sparse_input_weights
                    input_weights = sprand(population(pop_indx).nodes(i),  population(pop_indx).n_input_units+1, config.sparsity(i));
                    %input_weights(input_weights ~= 0) = ...
                    %    2*input_weights(input_weights ~= 0);
                else
                    input_weights = rand(population(pop_indx).nodes(i),  population(pop_indx).n_input_units+1);
                end
            case 'orth'
                input_weights = ones(population(pop_indx).nodes(i),  population(pop_indx).n_input_units+1);
        end
        population(pop_indx).input_weights{i} = input_weights;
        
        
        % add other necessary parameters
        % kernel
        population(pop_indx).pad_size(i) = randi([1 1]);
        population(pop_indx).kernel_stride(i) = config.vesicle_radius*2;%randi([1 5]);
        population(pop_indx).kernel_size(i) = config.vesicle_radius*2;%randi([1 5]);
        population(pop_indx).kernel{i} = ones(population(pop_indx).kernel_size(i))/population(pop_indx).kernel_size(i)^2; % summation filter
        
        % individual should keep track of final state for certain tasks
        population(pop_indx).last_state{i} = zeros(1,population(pop_indx).nodes(i));
    end
    
    %track total nodes in use
    population(pop_indx).total_units = 0;
    
    
    %% weights and connectivity of all reservoirs
    for i= 1:config.num_reservoirs
        
        for j= 1:config.num_reservoirs % only ensemble at the moment - no inter-reservoir connectivity
            
            % assign places for vesicles
            population(pop_indx).bitmatrix{i,j} = round(rand(population(pop_indx).vesicle_grid_size));
            
            % assign initial intensities
            phimatrix =zeros(size(population(pop_indx).bitmatrix{i,j}));
            phimatrix(population(pop_indx).bitmatrix{i,j} == 1) = rand(nnz(population(pop_indx).bitmatrix{i,j}),1)*0.06;
            population(pop_indx).phimatrix{i,j} = phimatrix;
            
            % source vesicles
            sourcematrix = population(pop_indx).bitmatrix{i,j};
            population(pop_indx).num_sources = randi([1 3]);%nnz(sourcematrix)]);
            % assign source locations
            possible_source_locations = zeros(nnz(sourcematrix),1);
            possible_source_locations(randperm(nnz(sourcematrix),population(pop_indx).num_sources)) = 1;
            sourcematrix(sourcematrix==1) = possible_source_locations;
            population(pop_indx).sourcematrix{i,j} = sourcematrix;     
        end
        % count total nodes including sub-reservoir nodes
        population(pop_indx).total_units = population(pop_indx).total_units + population(pop_indx).nodes(i);
    end
    
    
    % Add random output weights - these are typically trained for tasks but
    % can be evolved as well
    if config.add_input_states
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units + population(pop_indx).n_input_units, population(pop_indx).n_output_units)-1;
    else
        population(pop_indx).output_weights = 2*rand(population(pop_indx).total_units, population(pop_indx).n_output_units)-1;
    end
    
    % Add placeholder for behaviours
    population(pop_indx).behaviours = [];
    
end
