%% Mutation operator used for magnetic film

function offspring = mutateMultiMM(offspring,config)

for layer = 1:config.num_layers
    
    input_scaling = offspring.layer(layer).input_scaling(:);
    pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
    input_scaling(pos) = mutateWeight(input_scaling(pos),[-config.input_scaler, config.input_scaler],config);
    offspring.layer(layer).input_scaling = reshape(input_scaling,size(offspring.layer(layer).input_scaling));
    
    leak_rate = offspring.layer(layer).leak_rate(:);
    pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
    leak_rate(pos) = mutateWeight(leak_rate(pos),[0, 1],config);
    offspring.layer(layer).leak_rate = reshape(leak_rate,size(offspring.layer(layer).leak_rate));
    
    % vampire params
    damping = offspring.layer(layer).damping(:);
    pos = randperm(length(damping),sum(rand(length(damping),1) < config.mut_rate));
    damping(pos) = mutateWeight(damping(pos),config.damping_parameter,config);
    offspring.layer(layer).damping = reshape(damping,size(offspring.layer(layer).damping));
    
    anisotropy = offspring.layer(layer).anisotropy(:);
    pos = randperm(length(anisotropy),sum(rand(length(anisotropy),1) < config.mut_rate));
    anisotropy(pos) = mutateWeight(anisotropy(pos),config.anisotropy_parameter,config);
    offspring.layer(layer).anisotropy = reshape(anisotropy,size(offspring.layer(layer).anisotropy));
    
    temperature = offspring.layer(layer).temperature(:);
    pos = randperm(length(temperature),sum(rand(length(temperature),1) < config.mut_rate));
    temperature(pos) = mutateWeight(temperature(pos),config.temperature_parameter,config);
    offspring.layer(layer).temperature = reshape(temperature,size(offspring.layer(layer).temperature));
    
    exchange = offspring.layer(layer).exchange(:);
    pos = randperm(length(exchange),sum(rand(length(exchange),1) < config.mut_rate));
    exchange(pos) = mutateWeight(exchange(pos),config.exchange_parameter,config);
    offspring.layer(layer).exchange = reshape(exchange,size(offspring.layer(layer).exchange));
    
    magmoment = offspring.layer(layer).magmoment(:);
    pos = randperm(length(magmoment),sum(rand(length(magmoment),1) < config.mut_rate));
    magmoment(pos) = mutateWeight(magmoment(pos),config.magmoment_parameter,config);
    offspring.layer(layer).magmoment = reshape(magmoment,size(offspring.layer(layer).magmoment));
    
    applied_field_strength = offspring.layer(layer).applied_field_strength(:);
    pos = randperm(length(applied_field_strength),sum(rand(length(applied_field_strength),1) < config.mut_rate));
    applied_field_strength(pos) = mutateWeight(applied_field_strength(pos),config.applied_field_strength,config);
    offspring.layer(layer).applied_field_strength = reshape(applied_field_strength,size(offspring.layer(layer).applied_field_strength));
    
    thickness = offspring.layer(layer).thickness(:);
    pos = randperm(length(thickness),sum(rand(length(thickness),1) < config.mut_rate));
    thickness(pos) = round(mutateWeight(thickness(pos),[0 1],config)*10)/10;
    offspring.layer(layer).thickness = reshape(thickness,size(offspring.layer(layer).thickness));
    
    time_steps_increment = offspring.layer(layer).time_steps_increment(:);
    pos = randperm(length(time_steps_increment),sum(rand(length(time_steps_increment),1) < config.mut_rate));
    time_steps_increment(pos) = round(mutateWeight(time_steps_increment(pos),[config.time_steps_increment(1) config.time_steps_increment(2)],config));
    offspring.layer(layer).time_steps_increment = reshape(time_steps_increment,size(offspring.layer(layer).time_steps_increment));
    
    if config.evolve_material_density
        material_density = offspring.layer(layer).material_density(:);
        pos = randperm(length(material_density),sum(rand(length(material_density),1) < config.mut_rate));
        material_density(pos) = mutateWeight(material_density(pos),[0 1],config);
        offspring.layer(layer).material_density = reshape(material_density,size(offspring.layer(layer).material_density));
    end
    
    % geometry
    if config.evolve_geometry
        if config.evolve_poly
            poly_coord = offspring.layer(layer).poly_coord(:);
            pos = randperm(length(poly_coord),sum(rand(length(poly_coord),1) < config.mut_rate));
            poly_coord(pos) = mutateWeight(poly_coord(pos),[config.lb 1],config);
            offspring.layer(layer).poly_coord = reshape(poly_coord,size(offspring.layer(layer).poly_coord));
        else
            geo_width = offspring.layer(layer).geo_width(:);
            pos = randperm(length(geo_width),sum(rand(length(geo_width),1) < config.mut_rate));
            geo_width(pos) = mutateWeight(geo_width(pos),[config.lb 1],config);
            offspring.layer(layer).geo_width = reshape(geo_width,size(offspring.layer(layer).geo_width));
            
            geo_height = offspring.layer(layer).geo_height(:);
            pos = randperm(length(geo_height),sum(rand(length(geo_height),1) < config.mut_rate));
            geo_height(pos) = mutateWeight(geo_height(pos),[config.lb 1],config);
            offspring.layer(layer).geo_height = reshape(geo_height,size(offspring.layer(layer).geo_height));
        end
    end
    
    minimum_weight = 0.01; % anything less than this will be set to zero
    prob_2_del = 0.25*~config.single_input;% if sigle input, set to zero
    
    
    %% cycle through all sub-reservoirs
    for i = 1:config.num_res_in_layer(layer)
        
        % input weights
        input_weights = offspring.layer(layer).input_weights{i};
        
        for n = 1:size(input_weights,1)
            indices = find(input_weights(n,:)); % find inputs in use
            pos = randperm(length(indices),sum(rand(length(indices),1) < config.mut_rate)); %get weights to delete
            
            if config.single_input
                % make rectangle
                xy = createRect(offspring.layer(layer).geo_height,offspring.layer(layer).geo_width);
                
                % % reset input weights
                square_film_dimensions = sqrt(offspring.layer(layer).nodes(i));
                [xq,yq] = meshgrid(linspace(0,1,square_film_dimensions),linspace(0,1,square_film_dimensions));
                [in,on] = inpolygon(xq,yq,xy(:,1),xy(:,2));
                inputs_in_use = find(in | on);
                
                if ~ismember(indices(pos),inputs_in_use)
                    input_weights(n,indices(pos)) = 0; % remove weight and reassign
                    input_loc = randperm(length(inputs_in_use),1);
                    pos = inputs_in_use(input_loc); % only single input to change
                else
                    pos = indices(pos);
                end
            else
                if rand < prob_2_del
                    input_weights(n,indices(pos)) = 0; % delete weight
                end
                pos = randperm(size(input_weights,2),length(pos)); % find new positions to change
            end
            input_weights(n,pos) = mutateWeight(input_weights(n,pos),[-1 1],config); % mutate any random weights, chance to mutate existing and non-existing weight
        
        input_weights(n,input_weights(n,:)<minimum_weight & input_weights(n,:)~=0 & input_weights(n,:)>-minimum_weight) = 0;
        
        end
        
        
%         if nnz(input_weights) > offspring.n_input_units + 1
%             error('Error: More than one input.\n')
%         end
        
        offspring.layer(layer).input_weights{i} = input_weights; %reshape(input_weights,size(offspring.layer(layer).input_weights{i}));
        
        % width of inputs
        if config.input_widths
            input_widths = offspring.layer(layer).input_widths{i}(:);
            pos =  randperm(length(input_widths),sum(rand(length(input_widths),1) < config.mut_rate));
            input_widths(pos) = ceil(abs(randn(length(pos),1)));
            input_widths(input_widths > round(sqrt(offspring.layer(layer).nodes(i))/4)) = round(sqrt(offspring.layer(layer).nodes(i))/4);% cap at 1/4 size of space
            offspring.layer(layer).input_widths{i} = reshape(input_widths,size(offspring.layer(layer).input_widths{i}));
        end
        
%         if config.random_alloy(i) || config.core_shell(i)
%             interfacial_exchange = offspring.layer(layer).interfacial_exchange(i);
%             pos = randperm(length(interfacial_exchange),sum(rand(length(interfacial_exchange),1) < config.mut_rate));
%             interfacial_exchange(pos) = mutateWeight(interfacial_exchange(pos),config.exchange_parameter,config);
%             offspring.layer(layer).interfacial_exchange(i) = reshape(interfacial_exchange,size(offspring.layer(layer).interfacial_exchange(i)));
%         end
%         
%         if config.random_alloy(i)
%             alloy_fraction = offspring.layer(layer).alloy_fraction(i);
%             pos = randperm(length(alloy_fraction),sum(rand(length(alloy_fraction),1) < config.mut_rate));
%             alloy_fraction(pos) = mutateWeight(alloy_fraction(pos),[0 1],config);
%             offspring.layer(layer).alloy_fraction = reshape(alloy_fraction,size(offspring.layer(layer).alloy_fraction));
%         end
%         
%         if config.core_shell(i)
%             shell_size = offspring.layer(layer).shell_size(i,2);
%             pos = randperm(length(shell_size),sum(rand(length(shell_size),1) < config.mut_rate));
%             shell_size(pos) = mutateWeight(shell_size(pos),[0 1],config);
%             offspring.layer(layer).shell_size(i,2) = reshape(shell_size,size(offspring.layer(layer).shell_size(i,2)));
%         end
        
        periodic_boundary = offspring.layer(layer).periodic_boundary(i,logical(config.periodic_boundary));
        pos = randperm(length(periodic_boundary),sum(rand(length(periodic_boundary),1) < config.mut_rate));
        periodic_boundary(pos) = round(mutateWeight(periodic_boundary(pos),[0 1],config));
        offspring.layer(layer).periodic_boundary(i,logical(config.periodic_boundary)) = reshape(periodic_boundary,size(offspring.layer(layer).periodic_boundary(i,logical(config.periodic_boundary))));      
    end
    
    % mutate output weights
    if config.evolve_output_weights
        output_weights = offspring.output_weights(:);
        indices = find(output_weights); % find outputs in use
        pos = randperm(length(indices),sum(rand(length(indices),1) < config.mut_rate)); %get weights to delete
        if rand < prob_2_del
            output_weights(indices(pos)) = 0; % delete weight
        end
        new_pos = randperm(length(output_weights),length(pos)); % find new positions to change
        output_weights(new_pos) = mutateWeight(output_weights(new_pos),[-config.output_weight_scaler config.output_weight_scaler],config); % mutate any random weights, chance to mutate existing and non-existing weight
        output_weights(output_weights<minimum_weight & output_weights~=0 & output_weights>-minimum_weight) = 0;
        offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
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
