%% mutateBZ.m
% Used to mutate BZ spcific parameters

function offspring = mutateBZ(offspring,config)

% params - input scaling and leak rate
input_scaling = offspring.input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),[-1, 1],config);
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0, 1],config);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

a = offspring.a(:);
pos =  randperm(length(a),sum(rand(length(a),1) < config.mut_rate));
a(pos) = mutateWeight(a(pos),[0, 1],config);
offspring.a = reshape(a,size(offspring.a));

b = offspring.b(:);
pos =  randperm(length(b),sum(rand(length(b),1) < config.mut_rate));
b(pos) = mutateWeight(b(pos),[0, 1],config);
offspring.b = reshape(b,size(offspring.b));

c = offspring.c(:);
pos =  randperm(length(c),sum(rand(length(c),1) < config.mut_rate));
c(pos) = mutateWeight(c(pos),[0, 1],config);
offspring.c = reshape(c,size(offspring.c));

time_period = offspring.time_period(:);
pos =  randperm(length(time_period),sum(rand(length(time_period),1) < config.mut_rate));
time_period(pos) = round(mutateWeight(time_period(pos),[1 config.max_time_period],config));
offspring.time_period = reshape(time_period,size(offspring.time_period));

input_length = offspring.input_length(:);
pos =  randperm(length(input_length),sum(rand(length(input_length),1) < config.mut_rate));
input_length(pos) = round(mutateWeight(input_length(pos),[1 offspring.time_period],config));
offspring.input_length = reshape(input_length,size(offspring.input_length));

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    for r = 1:3
        % input weights
        input_weights = offspring.input_weights{i}(:);
        pos =  randperm(length(input_weights),ceil(config.mut_rate*length(input_weights)));
        input_weights(pos) = mutateWeight(input_weights(pos),[-1, 1],config);
        offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
        
        % width of inputs
        if config.input_widths
            input_widths = offspring.input_widths{i,r}(:);
            pos =  randperm(length(input_widths),sum(rand(length(input_widths),1) < config.mut_rate));
            input_widths(pos) = ceil(abs(randn(length(pos),1)));
            input_widths(input_widths > round(sqrt(offspring.nodes(i))/4)) = round(sqrt(offspring.nodes(i))/4);% cap at 1/4 size of space
            offspring.input_widths{i,r} = reshape(input_widths,size(offspring.input_widths{i,r}));
        end
    end
    
end

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos =  randperm(length(output_weights),ceil(config.mut_rate*length(output_weights)));
    output_weights(pos) =  mutateWeight(output_weights(pos),[-1,1],config);
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end

% mutate feedback weights
if config.evolve_feedback_weights
    % feedback scaling
    feedback_scaling = offspring.feedback_scaling(:);
    pos =  randperm(length(feedback_scaling),sum(rand(length(feedback_scaling),1) < config.mut_rate));
    feedback_scaling(pos) = mutateWeight(feedback_scaling(pos),[-1,1],config);
    offspring.feedback_scaling = reshape(feedback_scaling,size(offspring.feedback_scaling));
    
    feedback_weights = offspring.feedback_weights(:);
    pos = randperm(length(feedback_weights),sum(rand(length(feedback_weights),1) < config.mut_rate));
       
    feedback_weights(pos) = mutateWeight(feedback_weights(pos),[-1 1],config);
    offspring.feedback_weights = reshape(feedback_weights,size(offspring.feedback_weights));
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
