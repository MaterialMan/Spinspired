%% Mutation operator used for wave-based reservoir systems
% Details:
%
function offspring = mutateWave(offspring,config)

% params - input scaling and leak rate
input_scaling = offspring.input_scaling(:);
pos =  randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),[-1,1],config);
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos =  randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0,1],config);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

% wave parameters
time_period = offspring.time_period(:);
pos =  randperm(length(time_period),sum(rand(length(time_period),1) < config.mut_rate));
time_period(pos) = round(mutateWeight(time_period(pos),[1, config.max_time_period],config));
offspring.time_period = reshape(time_period,size(offspring.time_period));

wave_speed = offspring.wave_speed(:);
pos =  randperm(length(wave_speed),sum(rand(length(wave_speed),1) < config.mut_rate));
wave_speed(pos) = round(mutateWeight(wave_speed(pos),[1,config.max_wave_speed],config));
offspring.wave_speed = reshape(wave_speed,size(offspring.wave_speed));

damping_constant = offspring.damping_constant(:);
pos =  randperm(length(damping_constant),sum(rand(length(damping_constant),1) < config.mut_rate));
damping_constant(pos) = mutateWeight(damping_constant(pos),[0,10],config);
offspring.damping_constant = reshape(damping_constant,size(offspring.damping_constant));

input_delay = offspring.input_delay(:);
pos =  randperm(length(input_delay),sum(rand(length(input_delay),1) < config.mut_rate));
input_delay(pos) = floor(mutateWeight(input_delay(pos),[1,config.max_input_delay],config));
offspring.input_delay = reshape(input_delay,size(offspring.input_delay));


% boundary_conditions
pos = randi([1 config.num_reservoirs]);
offspring.boundary_conditions(pos,:) = config.boundary_conditions{randi([1 length(config.boundary_conditions)])};

%define minimum weight
minimum_weight = 0.01; % anything less than this will be set to zero
prob_2_del = config.prune_rate;

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs

    % input weights
    input_weights = offspring.input_weights(:);
    pos =  randperm(length(input_weights),ceil(config.mut_rate*length(input_weights)));
    input_weights(pos) = mutateWeight(input_weights(pos),[-1,1],config);
    indices = find(input_weights); % find outputs in use
    del_pos = randperm(length(indices),sum(rand(length(indices),1) < prob_2_del)); %get weights to delete
    input_weights(indices(del_pos)) = 0; % delete weight
    input_weights(input_weights<minimum_weight & input_weights~=0 & input_weights>-minimum_weight) = 0; % remove small weights    
    offspring.input_weights = reshape(input_weights,size(offspring.input_weights));
    
    if config.input_widths
        input_widths = offspring.input_widths{i}(:);
        pos =  randperm(length(input_widths),sum(rand(length(input_widths),1) < config.mut_rate));
        input_widths(pos) = ceil(abs(randn(length(pos),1)));
        input_widths(input_widths > round(sqrt(offspring.nodes(i))/4)) = round(sqrt(offspring.nodes(i))/4);% cap at 1/4 size of space
        offspring.input_widths{i} = reshape(input_widths,size(offspring.input_widths{i}));
    end
    
    % hidden weights
    for j = 1:config.num_reservoirs
        switch(config.wave_system)
            case 'fully-connected'
                if i~=j
                    W = offspring.W{i,j}(:);
                    % select weights to change
                    pos =  randperm(length(W),ceil(config.mut_rate*length(W)));
                    W(pos) = mutateWeight(W(pos),[-1,1],config);
                    offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
                end
            case 'pipeline'
                if j == i+1
                    W = offspring.W{i,j}(:);
                    % select weights to change
                    pos =  randperm(length(W),ceil(config.mut_rate*length(W)));
                    W(pos) = mutateWeight(W(pos),[-1,1],config);
                    offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
                end
        end
    end
    
    offspring.connectivity(i,j) = nnz(offspring.W{i,j})/offspring.total_units.^2;
    
end

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos =  randperm(length(output_weights),ceil(config.mut_rate*length(output_weights)));
    output_weights(pos) =  mutateWeight(output_weights(pos),[-1,1],config);
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end

end

function value = mutateWeight(value,range,config)

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
