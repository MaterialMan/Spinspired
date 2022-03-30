function[final_states,individual] = collectMinRoRMTSplusStates(individual,input_sequence,config,target_output,dataset)

% check number of function inputs; unless specified data is seen as train
% data
if nargin < 5
    training_dataset =1;
else
    training_dataset =0;
end

% start from last state
if size(input_sequence,1) == 1
    temp_input_sequence = zeros(size(individual.last_state,1),size(input_sequence,2));
    temp_input_sequence(end,:) = input_sequence(end,:);
    input_sequence = temp_input_sequence;
    
    states = individual.last_state;
    start_n = size(individual.last_state,1)+1;
    end_n = size(individual.last_state,1)+1;
    
    single_input_task = 1;
else
    states = zeros(size(input_sequence,1),sum(config.num_nodes));
    start_n = 1;
    end_n = size(input_sequence,1);
    single_input_task = 0;
end

% preassign weights rather than calculate them iteratively
noise = randn(1,individual.total_units)*config.noise_level;
input = ([input_sequence repmat(individual.bias_node,size(input_sequence,1),1)]*(individual.input_weights.*individual.input_scaling));
W = (individual.W.*individual.W_scaling);
if config.evolve_feedback_weights
    W_fb = individual.feedback_scaling*individual.feedback_weights';
end

if (config.evolve_feedback_weights || config.teacher_forcing) && exist('target_output','var')
    target_output = target_output + config.noise_ratio*rand(size(target_output)); % add noise to regulate weights
    force_input = ([target_output repmat(individual.bias_node,size(target_output,1),1)]*(individual.input_weights.*individual.input_scaling));
end

%equation: x(n) = f(Win*u(n) + S)
for n = start_n:end_n
    
    if n > config.max_update_cycle
        
        if config.evolve_feedback_weights || config.teacher_forcing
            y_out = states(n-1,:)*individual.output_weights;
        end
        
        % calculate delay states
         % get delayed states
        if config.leak_delay
            tmp_states = states(n-1,:);
        else
            tmp_states = sum(individual.W_delay.*states(n-config.max_update_cycle:n-1,:),1);
        end
                                        
        if config.multi_activ
            % calculate states with specific activation functions
            for p = 1:length(config.activ_list)  
                prev_state = tmp_states*W;
                indx = individual.activ_Fcn_indx == p;
                if sum(indx) > 0
                    states(n,indx) = config.activ_list{p}(input(n-1,indx) + prev_state(indx) + noise(indx));
                end
            end
        else
            %feeds trained output to the input (no feedback weights!)
            %feeds trained output to the input (no feedback weights!)
            if config.teacher_forcing
                if ~training_dataset % only apply generative mode to test data
                    if n < config.generate_n  % not in generative mode
                        if config.evolve_feedback_weights
                            states(n,:) = config.activ_list{1}(input(n-1,:) + tmp_states*W + noise + target_output(n,:)*W_fb);
                        else
                            states(n,:) = config.activ_list{1}(force_input(n-1,:) + tmp_states*W + noise);
                        end
                    else % in generative mode
                        if config.evolve_feedback_weights
                            states(n,:) = config.activ_list{1}(input(n-1,:) + tmp_states*W + noise + y_out*W_fb);
                        else
                            states(n,:) = config.activ_list{1}(([y_out individual.bias_node]*(individual.input_weights.*individual.input_scaling)) + tmp_states*W + noise);
                        end
                    end 
                else % collect data as normal for training 
                    states(n,:) = config.activ_list{1}(force_input(n-1,:) + tmp_states*W + noise);
                end
            else
                if config.evolve_feedback_weights %no training occurs, weights are evolved
                    states(n,:) = config.activ_list{1}(input(n-1,:) + tmp_states*W + noise + y_out*W_fb);
                else
                    % default state collection
                    states(n,:) = config.activ_list{1}(input(n-1,:) + tmp_states*W + noise);
                end
            end
        end
        
        % add leak states
        if config.leak_on && ~config.post_leak_on
            if config.leak_delay
                tmp_states = sum(individual.W_delay.*states(n-config.max_update_cycle:n-1,:));
                states(n,:) = (1-individual.leak_rate).*tmp_states + individual.leak_rate.*states(n,:);
            else
                states(n,:) = (1-individual.leak_rate).*states(n-1,:) + individual.leak_rate.*states(n,:);
            end
        end
        
    end
end

% add IIR filter
%     if config.iir_filter_on && n > config.iir_filter_order
%         prev_states = states{i}(n-(1:config.iir_filter_order),:);
%         filter_states(n+1,:) = iirFilterStates(prev_states,filter_states(n-(1:config.iir_filter_order),:),individual,config);
%     end


%assign last state variable
individual.last_state(1,:) = []; % pop
individual.last_state = [individual.last_state; states(end,:)];

% concat input states
if config.add_input_states == 1
    final_states = [states input_sequence];
else
    final_states = states;
end

if size(input_sequence,1) == config.max_update_cycle+2
    final_states = final_states(end,:); % remove washout
else
    final_states = final_states(config.wash_out+1:end,:); % remove washout
end

if single_input_task
    final_states = final_states(end,:);
end

