%% assess_ReservoirName_.m
% Template function to collect reservoir states. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
%
% This is called by the @config.assessFcn pointer.

function[final_states,individual] = assessHeterotic(individual,input_sequence,config,target_output)

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
    x{i} = zeros(size(input_sequence,1),individual.nodes(i));
end

% pre-assign anything that can be calculated before running the reservoir
%.....
switch(config.architecture)
    
    case 'RoR'
        
        % run through input sequence
        for n = 2:size(input_sequence,1)
            
            % Calculate reservoir states - general state equation for multi-reservoir system: x(n) = f(Win*u(n) + S)
            for i= 1:config.num_reservoirs
                
                % gather previous states of all sub reservoirs
                for k= 1:config.num_reservoirs
                    if i==k % other sub reservoirs
                        x{i}(n,:) = x{i}(n,:) + states{k}(n-1,:);
                    else
                        x{i}(n,:) = x{i}(n,:) + individual.W_switch(i,k)*((individual.W{i,k}*individual.W_scaling(i,k))*states{k}(n-1,:)')';
                    end
                end
                
                % adjust states of subresevoir to take into account other
                % subreservoirs
                individual.res{i}.last_state{1} = x{i}(n,:);
                
                input = input_sequence(n,:);
%                 if individual.subres_config{i}.add_input_states
%                     input =  [states{i-1} input_sequence];%* (individual.W_scaling(i-1,i)*individual.W{i-1,i});
%                 else
%                     input =  [states{i-1}];
%                 end
                
                individual.subres_config{i}.wash_out = 0;
                [state,individual.res{i}]= individual.subres_config{i}.assessFcn(individual.res{i},input,individual.subres_config{i});
                states{i}(n,:) = state;
            end
            
        end
        
    case 'ensemble'
        
        for i= 1:config.num_reservoirs
            individual.subres_config{i}.wash_out = 0;
            
            input = input_sequence;
            
            states{i}= individual.subres_config{i}.assessFcn(individual.res{i},input,individual.subres_config{i});
        end
        
    case {'pipeline','pipeline_IA'}
        
        for i= 1:config.num_reservoirs
            individual.subres_config{i}.wash_out = 0;
            
            if i < 2
                input = input_sequence;
            else
                if individual.subres_config{i}.add_input_states
                    input =  [states{i-1} input_sequence];%* (individual.W_scaling(i-1,i)*individual.W{i-1,i});
                else
                    input =  [states{i-1}];
                end
            end
            
            temp_indv = individual;
            temp_indv.subres_config{i}.add_input_states = 0;
            states{i}= individual.subres_config{i}.assessFcn(temp_indv.res{i},input,temp_indv.subres_config{i});
        end
        
end

% Add leak states, if used
if config.leak_on
    states= getLeakStates(states,individual,input_sequence,config);
end

% Concat all states for output weights
colors = distinguishable_colors(config.num_reservoirs);
final_states = [];
for i= 1:config.num_reservoirs
    final_states = [final_states states{i}];
    
    %assign last state variable
    individual.last_state{i} = states{i}(end,:);
    
    if config.plot_states
        plot(states{i},'Color',colors(i,:))
        hold on
    end
end
hold off
%legend
drawnow

%pause(0.5)

% Concat input states
if config.add_input_states == 1
    final_states = [final_states input_sequence];
end

% Remove washout and output final states
final_states = final_states(config.wash_out+1:end,:);

