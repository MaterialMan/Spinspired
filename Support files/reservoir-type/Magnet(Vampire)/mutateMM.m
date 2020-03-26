%% Mutation operator used for magnetic film

function offspring = mutateMM(offspring,config)

input_scaling = offspring.input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = input_scaling(pos) + randn*0.15;%2*rand(length(pos),1)-1;
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = leak_rate(pos) + randn*0.15; %rand(length(pos),1);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

% vampire params
damping = offspring.damping(:);
pos = randperm(length(damping),sum(rand(length(damping),1) < config.mut_rate));
damping(pos) = mutateWeight(damping(pos),config.damping_parameter,config);
offspring.damping = reshape(damping,size(offspring.damping));

anisotropy_parameter = offspring.anisotropy_parameter(:);
pos = randperm(length(anisotropy_parameter),sum(rand(length(anisotropy_parameter),1) < config.mut_rate));
anisotropy_parameter(pos) = mutateWeight(anisotropy_parameter(pos),config.anisotropy_parameter,config);
offspring.anisotropy_parameter = reshape(anisotropy_parameter,size(offspring.anisotropy_parameter));

temperature_parameter = offspring.temperature_parameter(:);
pos = randperm(length(temperature_parameter),sum(rand(length(temperature_parameter),1) < config.mut_rate));
temperature_parameter(pos) = mutateWeight(temperature_parameter(pos),config.temperature_parameter,config);
offspring.temperature_parameter = reshape(temperature_parameter,size(offspring.temperature_parameter));

exchange_parameter = offspring.exchange_parameter(:);
pos = randperm(length(exchange_parameter),sum(rand(length(exchange_parameter),1) < config.mut_rate));
exchange_parameter(pos) = mutateWeight(exchange_parameter(pos),config.exchange_parameter,config);
offspring.exchange_parameter = reshape(exchange_parameter,size(offspring.exchange_parameter));

magmoment_parameter = offspring.magmoment_parameter(:);
pos = randperm(length(magmoment_parameter),sum(rand(length(magmoment_parameter),1) < config.mut_rate));
magmoment_parameter(pos) = mutateWeight(magmoment_parameter(pos),config.magmoment_parameter,config);
offspring.magmoment_parameter = reshape(magmoment_parameter,size(offspring.magmoment_parameter));

applied_field_strength = offspring.applied_field_strength(:);
pos = randperm(length(applied_field_strength),sum(rand(length(applied_field_strength),1) < config.mut_rate));
applied_field_strength(pos) = mutateWeight(applied_field_strength(pos),config.applied_field_strength,config);
offspring.applied_field_strength = reshape(applied_field_strength,size(offspring.applied_field_strength));


%% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
    input_weights = offspring.input_weights{i}(:);
    pos =  randperm(length(input_weights),ceil(config.mut_rate*length(input_weights)));
    for n = 1:length(pos)
        input_weights(pos(n)) = input_weights(pos(n)) - randn*0.15;
    end
    offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
    
    % width of inputs
    if config.input_widths
        input_widths = offspring.input_widths{i}(:);
        pos =  randperm(length(input_widths),sum(rand(length(input_widths),1) < config.mut_rate));
        input_widths(pos) = ceil(abs(randn(length(pos),1)));
        input_widths(input_widths > round(sqrt(offspring.nodes(i))/4)) = round(sqrt(offspring.nodes(i))/4);% cap at 1/4 size of space
        offspring.input_widths{i} = reshape(input_widths,size(offspring.input_widths{i}));
    end
    
    % hidden weights
    for j = 1:config.num_reservoirs
        switch(offspring.architecture)
            case 'pipeline'
                if j == i+1
                    W = offspring.W{i,j}(:);
                    % select weights to change
                    pos =  randperm(length(W),ceil(config.mut_rate*length(W)));
                    for n = 1:length(pos)
                        W(pos(n)) = mutateWeight(W(pos(n)),[-1, 1],config);
                    end
                    offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
                end
        end
        
         offspring.connectivity(i,j) = nnz(offspring.W{i,j})/offspring.total_units.^2;
    end      
end

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos =  randperm(length(output_weights),ceil(config.mut_rate*length(output_weights)));
    
    for n = 1:length(pos)
      output_weights(pos(n)) =  mutateWeight(output_weights(pos(n)),[-1,1],config);
    end
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end

end

function value = mutateWeight(value,range,config)

switch(config.mutate_type)
    case 'gaussian'
        flag = 1;
        while(flag)
            t_value = value + (range(1) + (range(2)-range(1))*randn);
            
            % check within range
            if t_value <= range(2) && t_value >= range(1)
                flag = 0;
            end
        end
        value = t_value;
        
    case 'uniform'
        value = range(1) + (range(2)-range(1))*rand;
end
end
