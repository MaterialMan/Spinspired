%% assess_ReservoirName_.m
% Template function to collect reservoir states. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
%
% This is called by the @config.assessFcn pointer.

function[final_states,individual] = assessIsing(individual,input_sequence,config,target_output)

%if single input entry, add previous state
if size(input_sequence,1) == 1
    input_sequence = [zeros(size(input_sequence)); input_sequence];
end

% pre-allocate state matrices
for i= 1:config.num_reservoirs
    if size(input_sequence,1) == 2
        states{i} = individual.last_state{i};
    else
        states{i} = zeros(size(input_sequence,1),individual.nodes(i));
    end

    node_grid_size(i) = config.num_nodes(i);
    
    %% preassign allocate input sequence and time multiplexing
    input{i} = [input_sequence repmat(individual.bias_node(i),size(input_sequence,1),1)]*(individual.input_weights{i}*individual.input_scaling(i))';
    
    % time multiplex -
    input_mul{i} = zeros(size(input_sequence,1)*individual.time_period(i),size(input{i},2));
    if individual.time_period > 1
        input_mul{i}(mod(1:size(input_mul{i},1),individual.time_period(i)) == 1,:) = input{i};
    else
        input_mul{i} = input{i};
    end
    
    % change input widths
    for n = 1:size(input_mul{i},1)
        m = reshape(input_mul{i}(n,:),node_grid_size(i),node_grid_size(i));
        f_pos = find(m);
        input_matrix_2d = m;
        for p = 1:length(f_pos)
            t = zeros(size(m));
            t(f_pos(p)) = m(f_pos(p));
            [t] = adjustInputShape(t,individual.input_widths{i}(f_pos(p)));
            input_matrix_2d = input_matrix_2d + t;
        end
        input_mul{i}(n,:) = input_matrix_2d(:);
    end       
    
    % pre-assign anything that can be calculated before running the reservoir
    spin{i} = individual.init_spin{i};
end

% Calculate reservoir states - general state equation for multi-reservoir system: x(n) = f(Win*u(n) + S)
for n = 2:size(input_mul{i},1)
    
    for i= 1:config.num_reservoirs % cycle through sub-reservoirs
        
%         for k= 1:config.num_reservoirs % collect previous states of all sub-reservoirs and multiply them by the connecting matrices `W`
%             x{i}(n,:) = x{i}(n,:) + ((individual.W{i,k}*individual.W_scaling(i,k))*states{k}(n-1,:)')';
%         end

        % assign input
        spin{i}(input_mul{i}(n,:) ~= 0) = sign(input_mul{i}(n,(input_mul{i}(n,:) ~= 0)));
        
        %spin{i}(input_mul{i}(n,:) ~= 0) = -spin{i}(input_mul{i}(n,:) ~= 0);
        
        % calulate states
        spin{i} = metropolis(spin{i}, individual.kT, individual.J);
        states{i}(n,:) = reshape(spin{i},1,config.num_nodes(i).^2);
    end
end

% Add leak states, if used
if config.leak_on
    states = getLeakStates(states,individual,input_sequence,config);
end

% Concat all states for output weights
final_states = [];
for i= 1:config.num_reservoirs
    final_states = [final_states states{i}];
    
    %assign last state variable
    individual.last_state{i} = states{i}(end,:);
end

% Concat input states
if config.add_input_states == 1
    final_states = [final_states input_sequence];
end

% Remove washout and output final states
if size(input_sequence,1) == 2
    final_states = final_states(end,:); % remove washout
else
    final_states = final_states(config.wash_out+1:end,:); % remove washout
end