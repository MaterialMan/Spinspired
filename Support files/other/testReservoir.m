function [individual,test_states,test_sequence] = testReservoir(individual,config)

switch(config.input_mechanism)
    
    case 'BSA'
        %encode
        [encoded_data,~,filter,train_output_sequence]= encodeBSAspike(config.train_input_sequence, config.filter,config.train_output_sequence);
        %run
        wash_out = config.wash_out;
        config.wash_out = 0;
        train_states = config.assessFcn(individual,encoded_data,config,config.train_output_sequence);
        config.wash_out = wash_out;
        
        %decode
        train_states = decodeBSAspike(filter,(train_states),config.wash_out);
        
        
        train_output_sequence =config.train_output_sequence;
        
        if config.val_fraction > 0
            %encode
            [encoded_data,~,filter,val_output_sequence]= encodeBSAspike(config.val_input_sequence, config.filter,config.val_output_sequence);
            %run
            wash_out = config.wash_out;
            config.wash_out = 0;
            val_states = config.assessFcn(individual,encoded_data,config,config.val_output_sequence);
            config.wash_out = wash_out;
            
            %decode
            val_states = decodeBSAspike(filter,(val_states),config.wash_out);
            
            val_output_sequence = config.val_output_sequence;
        end
        
    otherwise
        train_output_sequence =config.train_output_sequence;
        train_states = config.assessFcn(individual,config.train_input_sequence,config,config.train_output_sequence);
        if config.val_fraction > 0
            val_output_sequence = config.val_output_sequence;
            val_states = config.assessFcn(individual,config.val_input_sequence,config,config.val_output_sequence);
        end
end

%if W_out are evolved instead of trained
if config.evolve_output_weights
    output_train_sequence = config.outputFcn(train_states*individual.output_weights);
    individual.train_error = calculateError(output_train_sequence,config.train_output_sequence,config);
    if config.val_fraction > 0
        output_val_sequence = config.outputFcn(val_states*individual.output_weights);
        individual.val_error = calculateError(output_val_sequence,config.val_output_sequence,config);
    else
        individual.val_error = 0;
    end
else
    % Find best reg parameter
    reg_train_error = [];
    reg_val_error =[];
    reg_weights=[];
    reg_param = [10e-1 10e-3 10e-5 10e-7 10e-9 10e-11];
    
    % assign output function, pre-training
    %train_states = config.outputFcn(train_states);
    
    for i = 1:length(reg_param)
        %Train: tanspose is inversed compared to equation
        output_weights = train_output_sequence(config.wash_out+1:end,:)'*train_states*inv(train_states'*train_states + reg_param(i)*eye(size(train_states'*train_states)));
        
        % Calculate trained output Y
        output_train_sequence = train_states*output_weights';
        reg_train_error(i,:)  = calculateError(output_train_sequence,train_output_sequence,config);
        
        % Calculate trained output Y
        if config.val_fraction > 0
            output_val_sequence = val_states*output_weights';
            reg_val_error(i,:)  = calculateError(output_val_sequence,val_output_sequence,config);
        end
        reg_weights(i,:,:) =output_weights';
    end
    
    if config.val_fraction > 0
        [~, reg_indx]= min(sum(reg_val_error,2));
        individual.val_error = sum(reg_val_error(reg_indx,:));
    else
        [~, reg_indx]= min(sum(reg_train_error,2));
        individual.val_error = 0;
    end
    individual.train_error = sum(reg_train_error(reg_indx,:));
        
    if config.identity_readout
        individual.output_weights = (eye(size(train_states,2)));
    else
        individual.output_weights = reshape(reg_weights(reg_indx,:,:),size(reg_weights,2),size(reg_weights,3));
    end
    
    %remove NaNs
    individual.output_weights(isnan(individual.output_weights)) = 0;
end

%% Evaluate on test data
if config.test_fraction > 0
    switch(config.input_mechanism)
        
        case 'spiking'
            %encode
            [encoded_data,~,filter,test_output_sequence]= encodeBSAspike(config.test_input_sequence, config.filter,config.test_output_sequence);
            %run
            wash_out = config.wash_out;
            config.wash_out = 0;
            test_states = config.assessFcn(individual,encoded_data,config,config.test_output_sequence,'test');
            config.wash_out = wash_out;
            %decode
            test_states = decodeBSAspike(filter,(test_states),config.wash_out);
            test_output_sequence =config.test_output_sequence;
            
        otherwise
            test_output_sequence =config.test_output_sequence;
            test_states = config.assessFcn(individual,config.test_input_sequence,config,test_output_sequence,'test');
    end
    
    test_sequence = config.outputFcn(test_states*individual.output_weights);
    individual.test_error = sum(calculateError(test_sequence,test_output_sequence,config,'test'));
else
    individual.test_error = 0;
end