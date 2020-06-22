%% Mutation operator used for magnetic film

function offspring = mutateMM(offspring,config)

input_scaling = offspring.input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),[-config.input_scaler, config.input_scaler],config);
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0, 1],config);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

% vampire params
damping = offspring.damping(:);
pos = randperm(length(damping),sum(rand(length(damping),1) < config.mut_rate));
damping(pos) = mutateWeight(damping(pos),config.damping_parameter,config);
offspring.damping = reshape(damping,size(offspring.damping));

anisotropy = offspring.anisotropy(:);
pos = randperm(length(anisotropy),sum(rand(length(anisotropy),1) < config.mut_rate));
anisotropy(pos) = mutateWeight(anisotropy(pos),config.anisotropy_parameter,config);
offspring.anisotropy = reshape(anisotropy,size(offspring.anisotropy));

temperature = offspring.temperature(:);
pos = randperm(length(temperature),sum(rand(length(temperature),1) < config.mut_rate));
temperature(pos) = mutateWeight(temperature(pos),config.temperature_parameter,config);
offspring.temperature = reshape(temperature,size(offspring.temperature));

exchange = offspring.exchange(:);
pos = randperm(length(exchange),sum(rand(length(exchange),1) < config.mut_rate));
exchange(pos) = mutateWeight(exchange(pos),config.exchange_parameter,config);
offspring.exchange = reshape(exchange,size(offspring.exchange));

magmoment = offspring.magmoment(:);
pos = randperm(length(magmoment),sum(rand(length(magmoment),1) < config.mut_rate));
magmoment(pos) = mutateWeight(magmoment(pos),config.magmoment_parameter,config);
offspring.magmoment = reshape(magmoment,size(offspring.magmoment));

applied_field_strength = offspring.applied_field_strength(:);
pos = randperm(length(applied_field_strength),sum(rand(length(applied_field_strength),1) < config.mut_rate));
applied_field_strength(pos) = mutateWeight(applied_field_strength(pos),config.applied_field_strength,config);
offspring.applied_field_strength = reshape(applied_field_strength,size(offspring.applied_field_strength));

if config.evolve_material_density
    material_density = offspring.material_density(:);
    pos = randperm(length(material_density),sum(rand(length(material_density),1) < config.mut_rate));
    material_density(pos) = mutateWeight(material_density(pos),[0 1],config);
    offspring.material_density = reshape(material_density,size(offspring.material_density));
end

thickness = offspring.thickness(:);
pos = randperm(length(thickness),sum(rand(length(thickness),1) < config.mut_rate));
thickness(pos) = round(mutateWeight(thickness(pos),[0 1],config)*10)/10;
offspring.thickness = reshape(thickness,size(offspring.thickness));

time_steps_increment = offspring.time_steps_increment(:);
pos = randperm(length(time_steps_increment),sum(rand(length(time_steps_increment),1) < config.mut_rate));
time_steps_increment(pos) = round(mutateWeight(time_steps_increment(pos),[config.time_steps_increment(1) config.time_steps_increment(2)],config));
offspring.time_steps_increment = reshape(time_steps_increment,size(offspring.time_steps_increment));

%% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
    input_weights = offspring.input_weights{i}(:);
    pos =  randperm(length(input_weights),ceil(config.mut_rate*length(input_weights)));
    input_weights(pos) = mutateWeight(input_weights(pos),[-config.input_scaler, config.input_scaler],config);
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
                    W(pos) = mutateWeight(W(pos),[-1, 1],config);
                    offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
                end
        end
        
        offspring.connectivity(i,j) = nnz(offspring.W{i,j})/offspring.total_units.^2;
    end
    
    
    if config.random_alloy(i) || config.core_shell(i)
        interfacial_exchange = offspring.interfacial_exchange(i);
        pos = randperm(length(interfacial_exchange),sum(rand(length(interfacial_exchange),1) < config.mut_rate));
        interfacial_exchange(pos) = mutateWeight(interfacial_exchange(pos),config.exchange_parameter,config);
        offspring.interfacial_exchange(i) = reshape(interfacial_exchange,size(offspring.interfacial_exchange(i)));
    end
    
    if config.random_alloy(i)
        alloy_fraction = offspring.alloy_fraction(i);
        pos = randperm(length(alloy_fraction),sum(rand(length(alloy_fraction),1) < config.mut_rate));
        alloy_fraction(pos) = mutateWeight(alloy_fraction(pos),[0 1],config);
        offspring.alloy_fraction = reshape(alloy_fraction,size(offspring.alloy_fraction));
    end
    
    if config.core_shell(i)
        shell_size = offspring.shell_size(i,2);
        pos = randperm(length(shell_size),sum(rand(length(shell_size),1) < config.mut_rate));
        shell_size(pos) = mutateWeight(shell_size(pos),[0 1],config);
        offspring.shell_size(i,2) = reshape(shell_size,size(offspring.shell_size(i,2)));
    end
    
    periodic_boundary = offspring.periodic_boundary(i,logical(config.periodic_boundary));
    pos = randperm(length(periodic_boundary),sum(rand(length(periodic_boundary),1) < config.mut_rate));
    periodic_boundary(pos) = round(mutateWeight(periodic_boundary(pos),[0 1],config));
    offspring.periodic_boundary(i,logical(config.periodic_boundary)) = reshape(periodic_boundary,size(offspring.periodic_boundary(i,logical(config.periodic_boundary))));
    
    
end

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos =  randperm(length(output_weights),ceil(config.mut_rate*length(output_weights)));
    output_weights(pos) =  mutateWeight(output_weights(pos),[-config.output_weight_scaler,config.output_weight_scaler],config);
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end

end

function value = mutateWeight(value,range,config)

if range(1)~=range(2)
    switch(config.mutate_type)
        case 'gaussian'
            for i = 1:length(value)
                flag = 1;
                while(flag)
                    t_value = value(i) + (range(1) + (range(2)-range(1))*randn);
                    
                    % check within range
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
