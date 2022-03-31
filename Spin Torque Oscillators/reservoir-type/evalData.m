function [states] = evalData(individual,config,dataset)
% only used in the case of taking the final value - for example
% sequence-to-classification
if iscell(config.(strcat(dataset,'_input_sequence')))
    if strcmp(config.res_type,'multiMM') || strcmp(config.res_type,'STO')
        % length to reset and ignore
        wash_out = 10;
        batch_size = 100;
        %store orignal sequence
        orig_input = config.(strcat(dataset,'_input_sequence'));

        % par settings        
        num_batches = round(length(config.(strcat(dataset,'_input_sequence')))/batch_size);
        state_grab = []; batch_input = []; states =[];

        % setup input streams for batches
        for batch_num = 1:num_batches
            batch_point = (batch_num-1)*batch_size + 1;
            seq_len = [];input_seq = []; cnt = 0;
            % concat input sequence into batches
            for b = 0:batch_size-1
                try % look at valid cells
                    seq_len(b+1) = size(orig_input{b+batch_point},2);
                    inputtemp = orig_input{b+batch_point};
                    if size(inputtemp,2) == config.task_num_inputs
                        inputtemp = inputtemp';
                    end
                    input_seq =  [input_seq zeros(size(inputtemp,1),wash_out) inputtemp];
                    cnt = cnt + 1;
                end
            end

            batch_input{batch_num} = input_seq';
            state_grab{batch_num} = cumsum(seq_len + wash_out); % points to grab states for evaluation
        end

        func = config.assessFcn;
        for batch_num = 1:num_batches
            tmp_individual(batch_num) = individual;
            %reset core indx
            tmp_individual(batch_num).core_indx = batch_num;
        end

        parfor batch_num = 1:num_batches            
            % evalute reservoir on batch
            tmp_states = func(tmp_individual(batch_num),batch_input{batch_num},config,[]); 
            states = [states; tmp_states(state_grab{batch_num},:)];
	    tmp_states = [];
        end
    else
        % Runs cells and collect last state
        parfor batch_num = 1:length(config.(strcat(dataset,'_input_sequence')))
            input = config.(strcat(dataset,'_input_sequence')){batch_num}';
            if size(input,2) ~= config.task_num_inputs
                input = input';
            end
            tmp_states = config.assessFcn(individual,input,config,[]);
            states(batch_num,:) = tmp_states(end,:);
	    tmp_states = [];
        end
    end
else
    states = config.assessFcn(individual,config.(strcat(dataset,'_input_sequence')),config,config.(strcat(dataset,'_output_sequence')));
end
end