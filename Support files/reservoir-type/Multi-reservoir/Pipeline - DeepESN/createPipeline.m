function population = createPipeline(config)


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
    
    % rand number of inner ESN's
    for i = 1:config.num_reservoirs
        
        %define num of units
        population(pop_indx).nodes(i) = config.num_nodes(i);
        
        % Scaling and leak rate
        population(pop_indx).input_scaling(i) = 2*rand-1; %increases nonlinearity
        population(pop_indx).leak_rate(i) = rand;
                
        %inputweights
        if i == 1
            sparsity = config.sparsity;
            additional_input_units = 0;
        else
            sparsity = config.connecting_sparsity;
            additional_input_units = population(pop_indx).nodes(i-1);
        end
        
        switch(config.input_weight_initialisation)
            case 'norm' % normal distribution
                if config.sparse_input_weights
                    input_weights = sprandn(population(pop_indx).nodes(i), additional_input_units + population(pop_indx).n_input_units +1, sparsity);
                else
                    input_weights = randn(population(pop_indx).nodes(i),  additional_input_units+population(pop_indx).n_input_units+1);
                end
            case 'uniform' % uniform dist between -1 and 1
                if config.sparse_input_weights
                    input_weights = sprand(population(pop_indx).nodes(i),  additional_input_units+population(pop_indx).n_input_units+1, sparsity);
                    input_weights(input_weights ~= 0) = ...
                        2*input_weights(input_weights ~= 0)  - 1;
                else
                    input_weights = 2*rand(population(pop_indx).nodes(i),  additional_input_units+population(pop_indx).n_input_units+1)-1;
                end
            case 'orth'
                input_weights = ones(population(pop_indx).nodes(i),  additional_input_units+population(pop_indx).n_input_units+1);
        end
  
        population(pop_indx).input_weights{i} = input_weights;
        
        
        %assign different activations, if necessary
        if config.multi_activ
            activ_positions = randi(length(config.activ_list),1,population(pop_indx).nodes(i));
            for act = 1:length(activ_positions)
                population(pop_indx).activ_Fcn{i,act} = config.activ_list{activ_positions(act)};
            end
        else
            population(pop_indx).activ_Fcn = config.activ_list;
        end
        
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
            
            if i==j
                switch(config.internal_weight_initialisation)
                case 'norm' % normal distribution
                    internal_weights = sprandn(population(pop_indx).nodes(i), population(pop_indx).nodes(j), population(pop_indx).connectivity(i,j));
                   
                case 'uniform' % uniform dist between -1 and 1
                    internal_weights = sprand(population(pop_indx).nodes(i), population(pop_indx).nodes(j), population(pop_indx).connectivity(i,j));
                    internal_weights(internal_weights ~= 0) = ...
                        2*internal_weights(internal_weights ~= 0)  - 1;   
                case 'orth'
                    internal_weights = orth(rand(population(pop_indx).nodes(i), population(pop_indx).nodes(j)));
                case 'sparse_orth'
                    internal_weights = orth(full(sprand(population(pop_indx).nodes(i), population(pop_indx).nodes(j), population(pop_indx).connectivity(i,j))));   
                end
            
%                 internal_weights = sprand(population(pop_indx).nodes(i), population(pop_indx).nodes(i), population(pop_indx).connectivity(i,j));
%                 
%                 internal_weights(internal_weights ~= 0) = ...
%                     2*internal_weights(internal_weights ~= 0)  - 1;
%                 
%                 % assign scaling for inner weights
                population(pop_indx).W_scaling(i,j) = 2*rand-1;
                population(pop_indx).W{i,j} = internal_weights;
            else
                population(pop_indx).W{i,j} = 0;
            end
            
        end
        population(pop_indx).total_units = population(pop_indx).total_units + population(pop_indx).nodes(i);
    end
    
    
    % add rand output weights
    if config.add_input_states
        output_units = population(pop_indx).total_units + population(pop_indx).n_input_units;
    else
        output_units =population(pop_indx).total_units;
    end
    
    
    switch(config.internal_weight_initialisation)
        case 'norm' % normal distribution
            output_weights = config.output_weight_scaler.*sprandn(output_units, population(pop_indx).n_output_units, config.output_connectivity);
            
        case 'uniform' % uniform dist between -1 and 1
            output_weights = sprand(output_units, population(pop_indx).n_output_units, config.output_connectivity);
            output_weights(output_weights ~= 0) = ...
                2*output_weights(output_weights ~= 0)  - 0.5;
            output_weights = config.output_weight_scaler.*output_weights;
        case 'orth'
            output_weights = orth(rand(output_units, population(pop_indx).n_output_units));
        case 'sparse_orth'
            output_weights = orth(full(sprand(output_units, population(pop_indx).n_output_units, config.output_connectivity)));
    end
    population(pop_indx).output_weights = output_weights;
    
    % add rand feedback weights
    if config.evolve_feedback_weights
        population(pop_indx).feedback_scaling = 2*rand;
        switch(config.feedback_weight_initialisation)
            case 'norm' % normal distribution
                feedback_weights = sprandn(population(pop_indx).total_units, population(pop_indx).n_input_units, config.feedback_connectivity);
                
            case 'uniform' % uniform dist between -1 and 1
                feedback_weights = sprand(population(pop_indx).total_units, population(pop_indx).n_input_units, config.feedback_connectivity);
                feedback_weights(feedback_weights ~= 0) = ...
                    2*feedback_weights(feedback_weights ~= 0)  - 1;
            case 'orth'
                feedback_weights = orth(rand(population(pop_indx).total_units, population(pop_indx).n_input_units));
            case 'sparse_orth'
                feedback_weights = orth(full(sprand(population(pop_indx).total_units, population(pop_indx).n_input_units, config.feedback_connectivity)));
                
        end
        population(pop_indx).feedback_weights = feedback_weights;
    end
    
    population(pop_indx).behaviours = [];
end