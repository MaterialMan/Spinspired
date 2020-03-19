%% create_ReservoirName_.m
% Template function to define reservoir parameters. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
% Tip: Try maintain ordering as structs keep this ordering.

% This is called by the @config.createFcn pointer.

function population = createHeterotic(config)

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
    
    % iterate through subreservoirs
    for i = 1:config.num_reservoirs
        
        %define num of units
        population(pop_indx).nodes(i) = config.num_nodes(i);
        
        population(pop_indx).leak_rate(i) = rand;
        
        %assign different reservoirs
        % reset values
        config.figure_array = [];
        population(pop_indx).subres_config{i} = config; % copy all params
        population(pop_indx).subres_config{i}.res_type = config.res_type{i};
        population(pop_indx).subres_config{i}.num_nodes = population(pop_indx).nodes(i);
        population(pop_indx).subres_config{i}.pop_size = 1; % only create single sub-reservoir
        
        % get pointers
        population(pop_indx).subres_config{i} = selectReservoirType(population(pop_indx).subres_config{i});
        % get all res parameters
        population(pop_indx).subres_config{i} = getAdditionalParameters(population(pop_indx).subres_config{i});
        population(pop_indx).subres_config{i}.add_input_states = 0; % do not add inputs
        
        
        % create initial population
        population(pop_indx).res{i} = population(pop_indx).subres_config{i}.createFcn(population(pop_indx).subres_config{i});
        
        population(pop_indx).last_state{i} = zeros(1,population(pop_indx).nodes(i));       
        
        % individual should keep track of final state for certain tasks
        population(pop_indx).last_state{i} = zeros(1,population(pop_indx).nodes(i));
    end
    
    %track total nodes in use
    population(pop_indx).total_units = 0;
    
    
    %% weights and connectivity of all reservoirs
    for i= 1:config.num_reservoirs
        
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
                    %population(pop_indx).W{i,j} = zeros(population(pop_indx).nodes(i), population(pop_indx).n_input_units);
                    population(pop_indx).W{i,j} = zeros(population(pop_indx).res{i}.nodes, population(pop_indx).res{j}.nodes);
                    
                case {'pipeline','pipeline_IA'}  % pipeline structure
                    if  j == i+1
                        switch(config.internal_weight_initialisation)
                            case 'norm' % normal distribution
                                %internal_weights = sprandn(population(pop_indx).res{i}.nodes + population(pop_indx).res{j}.n_input_units, population(pop_indx).res{j}.n_input_units, population(pop_indx).connectivity(i,j));
                                internal_weights = sprandn(population(pop_indx).nodes(j),population(pop_indx).nodes(i)+ population(pop_indx).res{j}.n_input_units +1, population(pop_indx).connectivity(i,j));

                            case 'uniform' % uniform dist between -1 and 1
                                %internal_weights = sprand(population(pop_indx).res{i}.nodes + population(pop_indx).res{j}.n_input_units, population(pop_indx).res{j}.n_input_units, population(pop_indx).connectivity(i,j));
                                internal_weights = sprand(population(pop_indx).nodes(j),population(pop_indx).nodes(i)+ population(pop_indx).res{j}.n_input_units +1, population(pop_indx).connectivity(i,j));

                                internal_weights(internal_weights ~= 0) = ...
                                    2*internal_weights(internal_weights ~= 0)  - 1;
                        end
                        
                        % assign scaling for inner weights
                        population(pop_indx).W_scaling(i,j) = 4*rand-2;
                        population(pop_indx).W{i,j} = internal_weights;
                        population(pop_indx).res{j}.input_weights{1} = internal_weights;
                    else
                        % if self, do nothing
                        population(pop_indx).connectivity(i,j) = 0;
                        population(pop_indx).W_scaling(i,j) = 0;
                        %population(pop_indx).W{i,j} = zeros(population(pop_indx).nodes(i), population(pop_indx).n_input_units);
                        population(pop_indx).W{i,j} = zeros(population(pop_indx).res{i}.nodes, population(pop_indx).res{j}.nodes);
                    end
                    
                case {'RoR','RoR_IA'}
                    
                    population(pop_indx).subres_config{i}.leak_on = 0;
                    
                    if  i ~= j
                        switch(config.internal_weight_initialisation)
                            case 'norm' % normal distribution
                                %internal_weights = sprandn(population(pop_indx).nodes(i), population(pop_indx).n_input_units, population(pop_indx).connectivity(i,j));
                                internal_weights = sprandn(population(pop_indx).res{i}.nodes, population(pop_indx).res{j}.nodes, population(pop_indx).connectivity(i,j));
                                
                            case 'uniform' % uniform dist between -1 and 1
                                %internal_weights = sprand(population(pop_indx).nodes(i), population(pop_indx).n_input_units, population(pop_indx).connectivity(i,j));
                                internal_weights = sprand(population(pop_indx).res{i}.nodes, population(pop_indx).res{j}.nodes, population(pop_indx).connectivity(i,j));
                                
                                internal_weights(internal_weights ~= 0) = ...
                                    2*internal_weights(internal_weights ~= 0)  - 1;
                        end
                        % assign scaling for inner weights
                        population(pop_indx).W_scaling(i,j) = 4*rand-2;
                        population(pop_indx).W{i,j} = internal_weights;
                        
                    else
                        population(pop_indx).connectivity(i,j) = 1;
                        population(pop_indx).W_scaling(i,j) = 1;
                        %population(pop_indx).W{i,j} = zeros(population(pop_indx).nodes(i), population(pop_indx).n_input_units);
                        population(pop_indx).W{i,j} = ones(population(pop_indx).res{i}.nodes, population(pop_indx).res{j}.nodes);
                    end
                  
            end
        end
        
        % count total nodes including sub-reservoir nodes
        population(pop_indx).total_units = population(pop_indx).total_units + population(pop_indx).res{i}.nodes;
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

