function [individual,test_states,test_sequence] = evalReservoir(individual,config)

dataset_list = {'train','val','test'};
dataset_to_eval = [config.train_fraction > 0; config.val_fraction > 0; config.test_fraction > 0];

% magnetic systems
if strcmp(config.res_type,'multiMM') || strcmp(config.res_type,'STO') %%&& isempty(config.architecture)
     mag_system = 1;
else
    mag_system = 0;
end

% eval all datasets
for d = 1:length(dataset_list)
    test_individual(d) = individual;
    if mag_system
        test_individual(d).core_indx = d;
    end
end

% parallel data split evaluations
if (~iscell(config.train_input_sequence) && mag_system && ~config.plot_states) || config.task_parallel % not using batches, so eval datasets in parallel
     parfor d = 1:length(dataset_list)
        states{d} = [];
        if dataset_to_eval(d)
            states{d} = evalData(test_individual(d),config,dataset_list{d});
        end
    end
else %uses batched data - can use parafor in evalData
    for d = 1:length(dataset_list)
        states{d} = [];
        if dataset_to_eval(d)
            states{d} = evalData(test_individual(d),config,dataset_list{d});
        end
    end
end

% decode spikes
switch(config.input_mechanism)
    case 'spiking'
        %decode
        states{1} = decodeBSAspike(config.filter,(states{1}));
        states{2} = decodeBSAspike(config.filter,(states{2}));
        states{3} = decodeBSAspike(config.filter,(states{3}));
end

% get output weights
individual = train(individual,states{1},states{2},config);

%% Evaluate on test data
if config.test_fraction > 0
    test_sequence = states{3}*individual.output_weights;
    individual.test_error = sum(calculateError(test_sequence,config.test_output_sequence,config));
    test_states = states{3};
else
    individual.test_error = 0;
end   

end


function individual = train(individual,train_states,val_states,config)

% check if evolving weights
if config.evolve_output_weights
    config.training_type = 'evolve';
end

switch(config.training_type)

    case 'pinv'

        individual.output_weights = pinv(train_states)*config.train_output_sequence(config.wash_out+1:end,:);

    case 'ridge'
        % Find best reg parameter
        reg_train_error = [];
        reg_val_error =[];
        reg_weights=[];
        reg_param = [10e-1 10e-3 10e-5 10e-7 10e-9 10e-11];

        for i = 1:length(reg_param)
            %Train: tanspose is inversed compared to equation
            output_weights = config.train_output_sequence(config.wash_out+1:end,:)'*train_states*inv(train_states'*train_states + reg_param(i)*eye(size(train_states'*train_states)));

            % Calculate trained output Y
            output_train_sequence = train_states*output_weights';
            reg_train_error(i,:)  = calculateError(output_train_sequence,config.train_output_sequence,config);

            % Calculate trained output Y
            if config.val_fraction > 0
                output_val_sequence = val_states*output_weights';
                reg_val_error(i,:)  = calculateError(output_val_sequence,config.val_output_sequence,config);
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
        individual.output_weights =reshape(reg_weights(reg_indx,:,:),size(reg_weights,2),size(reg_weights,3));

        %remove NaNs
        individual.output_weights(isnan(individual.output_weights)) = 0;

    case 'evolve'
        % left blank
end

%calculate errors
output_train_sequence = train_states*individual.output_weights;
individual.train_error = calculateError(output_train_sequence,config.train_output_sequence,config);
if config.val_fraction > 0
    output_val_sequence = val_states*individual.output_weights;
    individual.val_error = calculateError(output_val_sequence,config.val_output_sequence,config);
end

end


