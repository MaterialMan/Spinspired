%% Mutation operator used for magnetic film

function offspring = mutateSTO(offspring,config)

for layer = 1:config.num_layers
    
% reservoir properties
offspring.layer(layer) = mutFcn(offspring.layer(layer),[-config.input_scaler, config.input_scaler],0,[],config,'input_scaling');
offspring.layer(layer) = mutFcn(offspring.layer(layer),[0, 1],0,[],config,'leak_rate');

%if config.hetero_exchange_interaction
offspring.layer(layer) = mutFcn(offspring.layer(layer),[-1, 1],0,[],config,'W_scaling');
%end

% material properties
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.damping_parameter,0,[],config,'damping');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.anisotropy_parameter,0,[],config,'anisotropy');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.temperature_parameter,0,[],config,'temperature');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.exchange_parameter,0,[],config,'exchange');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.magmoment_parameter,0,[],config,'magmoment');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.applied_field_strength,0,[],config,'applied_field_strength');

% timescales
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.interpolation_length,1,[],config,'interpolation_length');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.time_steps_increment,1,[],config,'time_steps_increment');

% spin parameters
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.spin_transfer_relaxation_torque,0,[],config,'spin_transfer_relaxation_torque');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.spin_transfer_precession_torque,0,[],config,'spin_transfer_precession_torque');
offspring.layer(layer) = mutFcn(offspring.layer(layer),config.spin_transfer_torque_asymmetry,0,[],config,'spin_transfer_torque_asymmetry');

%% cycle through all sub-reservoirs
for res_indx = 1:length(offspring.layer(layer).nodes)
    
    offspring.layer(layer) = mutFcn(offspring.layer(layer),[-1, 1],0,res_indx,config,'input_weights');
    
    if config.hetero_exchange_interaction
        offspring.layer(layer) = mutFcn(offspring.layer(layer),[-1, 1],0,res_indx,config,'input_weights');        
        % remove weights  
        offspring.layer(layer).W{res_indx}(~offspring.layer(layer).array_mask{res_indx}) = 0;
    end
end

% mutate output weights
if config.evolve_output_weights
    offspring = mutFcn(offspring,[-config.output_weight_scaler,config.output_weight_scaler],0,[],config,'output_weights');
end

end
end

function value = mutateWeight(value,range,config)

if range(1)~=range(2)
    switch(config.mutate_type)
        case 'gaussian'
            for i = 1:length(value)
                flag = 1;
                while(flag)
                    t_value = value(i) + (range(1) - (range(2)-range(1))*randn);
                    
                    % check within range
%                     if (abs(t_value) <= abs(range(2))) && (abs(t_value) >= abs(range(1)))
%                          flag = 0;
%                     end
                    
                    if (t_value <= range(2)) && (t_value >= range(1))
                        flag = 0;
                    end
                end
                value(i) = t_value;
            end
        case 'uniform'
            value = range(1) + (range(2)-range(1))*rand;
    end
end
end

function offspring = mutFcn(offspring,range,rnd,indx,config,parameter)
    
    if iscell(offspring.(parameter))
        tmp = offspring.(parameter){indx}(:);
    else
        tmp = offspring.(parameter)(:);
    end

    pos = randperm(length(tmp),sum(rand(length(tmp),1) < config.mut_rate));
    if rnd
        tmp(pos) = round(mutateWeight(tmp(pos),range,config));
    else
        tmp(pos) = mutateWeight(tmp(pos),range,config);
    end
    
    % delete if weight
    if contains(parameter,'weight')
        indices = find(tmp); % find outputs in use
        del_pos = randperm(length(indices),sum(rand(length(indices),1) < config.prune_rate)); %get weights to delete
        tmp(indices(del_pos)) = 0; % delete weight
    end
    
    if iscell(offspring.(parameter))
         offspring.(parameter){indx} = reshape(tmp,size(offspring.(parameter){indx}));
    else
         offspring.(parameter) = reshape(tmp,size(offspring.(parameter)));
    end
   
end