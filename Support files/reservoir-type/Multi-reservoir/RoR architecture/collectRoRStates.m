function[final_states,individual] = collectRoRStates(individual,input_sequence,config,target_output)

%if single input entry, add previous state
if size(input_sequence,1) == 1 % pad input
    input_sequence = [zeros(size(input_sequence)); input_sequence];
end

%states = cell; x = cell;
for i= 1:config.num_reservoirs
    if size(input_sequence,1) == 2
        states{i} = individual.last_state{i};
    else
        states{i} = zeros(size(input_sequence,1),individual.nodes(i));
    end
    x{i} = zeros(size(input_sequence,1),individual.nodes(i));
end

% preassign activation function calls
%index = cell;
if config.multi_activ
    for i= 1:config.num_reservoirs
        for p = 1:length(config.activ_list)
            index{i,p} = findActiv(individual.activ_Fcn(i,:),config.activ_list{p});
        end
    end
end

%equation: x(n) = f(Win*u(n) + S)
for n = 2:size(input_sequence,1)-1
    
    for i= 1:config.num_reservoirs
        
        % this part controls the effect of different reservoirs
        for k= 1:config.num_reservoirs
            % x{i}(n,:) = x{i}(n,:) + individual.W_switch(i,k)*((individual.W{i,k}*individual.W_scaling(i,k))*states{k}(n-1,:)')';  % (n-1)
            x{i}(n,:) = x{i}(n,:) + individual.W_switch(i,k)*((individual.W{i,k}*individual.W_scaling(i,k))*states{k}(n-1,:)')'; % (n)turn on/off sub-reservoirs, i.e. adds/removes forwards/feedback connectivity between sub-reservoirs
        end
        
        if config.multi_activ
            for p = 1:length(config.activ_list)
                if config.evolve_feedback_weights
                    states{i}(n,index{i,p}) = config.activ_list{p}(((individual.input_weights{i}(index{i,p},:)*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])')+ x{i}(n,index{i,p})' + individual.feedback_weights*states{i}(n-1,:)*individual.output_weights(1:end-config.add_input_states,:));
                else
                    states{i}(n,index{i,p}) = config.activ_list{p}(((individual.input_weights{i}(index{i,p},:)*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])')+ x{i}(n,index{i,p})');
                end
            end
        else
            if config.teacher_forcing %feeds trained output to the input (no feedback weights!)
                if sum(sum(input_sequence(n-1:n,:))) ~= 0 %input_sequence(n,:) ~= 0 && n <= size(input_sequence,1)/2%sum(sum(input_sequence(n-1:n,:))) ~= 0 % teacher forcing
                    if config.evolve_feedback_weights
                        target_output(n-1,:) = target_output(n-1,:) + config.noise_ratio*rand(size(target_output(n-1,:))); % add noise to regulate weights
                        input = ((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])');
                        y_out = target_output(n-1,:);%*individual.output_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:);
                        W_fb = individual.feedback_scaling*individual.feedback_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:);
                        %states{i}(n,:) = individual.activ_Fcn{i}(input + x{i}(n,:)' + W_fb*y_out');
                        states{i}(n,:) = individual.activ_Fcn{i}(input + x{i}(n,:)' + W_fb*y_out');
                    else
                        target_output(n-1,:) = target_output(n-1,:) + config.noise_ratio*rand(size(target_output(n-1,:))); % add noise to regulate weights
                        states{i}(n,:) = individual.activ_Fcn{i}(((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node target_output(n-1,:)])') + x{i}(n,:)');
                    end
                else
                    %  states{i}(n,:) = individual.activ_Fcn{i}(((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])') + x{i}(n,:)' + ((individual.feedback_scaling*individual.feedback_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:))*(states{i}(n-1,:)*individual.output_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:))'));
                    if config.evolve_feedback_weights
                        y_out = states{i}(n-1,:)*individual.output_weights; %(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:);
                        input = ((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])');
                        W_fb = individual.feedback_scaling*individual.feedback_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:);
                        states{i}(n,:) = individual.activ_Fcn{i}(input + x{i}(n,:)' + W_fb*y_out');
                    else
                        input = (states{i}(n-1,:)*individual.output_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:))';
                        states{i}(n,:) = individual.activ_Fcn{i}(((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node input])') + x{i}(n,:)');
                    end
                end
            else
                if config.evolve_feedback_weights %no training occurs, weights are evolved
                    input = (individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])';
                    y_out = states{i}(n-1,:)*individual.output_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:);
                    W_fb = individual.feedback_scaling*individual.feedback_weights(sum(individual.nodes(1:i-1))+1:sum(individual.nodes(1:i)),:);
                    states{i}(n,:) = individual.activ_Fcn{i}(input + x{i}(n,:)' + W_fb*y_out');
                else % normal operating reservoir - if RoR, first reservoir is the only one to take an input
                    %if i ==1 || config.add_input_states
                        states{i}(n,:) = individual.activ_Fcn{i}(((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])')+ x{i}(n,:)');
                        % states{i}(n,:) = individual.activ_Fcn{i}(((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node input_sequence(n,:)])')+ x{i}(n,:)');
                 
%                     else
%                         states{i}(n,:) = individual.activ_Fcn{i}(((individual.input_weights{i}*individual.input_scaling(i))*([individual.bias_node zeros(size(input_sequence(n,:)))])')+ x{i}(n,:)');
%                     end
                end
            end
            
        end
        
        if config.leak_on 
            states{i}(n+1,:) = (1-individual.leak_rate(i))*states{i}(n,:) + individual.leak_rate(i)*states{i}(n+1,:);
        end
    end
    
    % add IIR filter
%     if config.iir_filter_on && n > config.iir_filter_order
%         prev_states = states{i}(n-(1:config.iir_filter_order),:);
%         filter_states(n+1,:) = iirFilterStates(prev_states,filter_states(n-(1:config.iir_filter_order),:),individual,config);
%     end
    
    
end


% get leak states after collection
% if config.leak_on
%     states = getLeakStates(states,individual,input_sequence,config);
% end

% concat all states for output weights
final_states = [];
for i= 1:config.num_reservoirs
    final_states = [final_states states{i}];
    
    %assign last state variable
    individual.last_state{i} = states{i}(end,:);
end


% concat input states
if config.add_input_states == 1
    final_states = [final_states input_sequence];
end

if size(input_sequence,1) == 2
    final_states = final_states(end,:); % remove washout
else
    final_states = final_states(config.wash_out+1:end,:); % remove washout
end

