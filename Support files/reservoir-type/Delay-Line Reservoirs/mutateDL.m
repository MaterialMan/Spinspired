function offspring = mutateDL(offspring,config)

% params - input scaling and leak rate
input_scaling = offspring.input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),[-1, 1],config);
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0, 1],config);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

% DL  parameters
eta = offspring.eta(:);
pos = randperm(length(eta),sum(rand(length(eta),1) < config.mut_rate));
eta(pos) = mutateWeight(eta(pos),[0, 1],config);
offspring.eta = reshape(eta,size(offspring.eta));

gamma = offspring.gamma(:);
pos = randperm(length(gamma),sum(rand(length(gamma),1) < config.mut_rate));
gamma(pos) = mutateWeight(gamma(pos),[0, 1],config);
offspring.gamma = reshape(gamma,size(offspring.gamma));

p = offspring.p(:);
pos = randperm(length(p),sum(rand(length(p),1) < config.mut_rate));
p(pos) = round(mutateWeight(p(pos),[1, 20],config));
offspring.p = reshape(p,size(offspring.p));

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
     input_weights = offspring.input_weights{i}(:);
        pos =  randperm(length(input_weights),ceil(config.mut_rate*length(input_weights)));
        input_weights(pos) = mutateWeight(input_weights(pos),[-1, 1],config);
        offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));     
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
