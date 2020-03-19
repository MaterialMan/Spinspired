%% mutate_ReservoirName_.m
% Template function to mutate the offspring reservoir. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
%
% This is called by the @config.mutateFcn pointer.
%
% Additional details:
% - Number of weights mutated is based on mut_rate; 50% chance to change existing weight or remove it

function offspring = mutate_ReservoirName_(offspring,config)
     
% params - input scaling and leak rate
input_scaling = offspring.input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),config);
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),config);
offspring.input_scaling = reshape(leak_rate,size(offspring.leak_rate));

% mutate other parameters - template below
% temp_variable_name = offspring.parameter(:);
% pos =  randperm(length(temp_variable_name),sum(rand(length(temp_variable_name),1) < config.mut_rate));
% temp_variable_name (pos) = 2*rand(length(pos),1);
% offspring.parameter = reshape(temp_variable_name,size(offspring.parameter));


% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs

    % input weights
    input_weights = offspring.input_weights{i}(:);
    pos = randperm(length(input_weights),sum(rand(length(input_weights),1) < config.mut_rate));
    input_weights(pos) = mutateWeight(input_weights(pos),config); 
    offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
        
    % Add additional sub-reservoir specific changes
    % e.g., connection matrix 'W'
    %     for j = 1:config.num_reservoirs
    %         W = offspring.W{i,j}(:);
    %         % select weights to change
    %         pos =  randperm(length(W),ceil(config.mut_rate*length(W)));
    %         for n = 1:length(pos)
    %             if rand < 0.5 % 50% chance to zero weight
    %                 W(pos(n)) = 0;
    %             else
    %                 W(pos(n)) = 2*rand-1;
    %             end   
    %         end
    %         offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
    %     end
    
end

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos = randperm(length(output_weights),sum(rand(length(output_weights),1) < config.mut_rate));   
    output_weights(pos) = mutateWeight(output_weights(pos),config);
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end


% mutate feedback weights
if config.evolve_feedback_weights
    % feedback scaling
    feedback_scaling = offspring.feedback_scaling(:);
    pos =  randperm(length(feedback_scaling),sum(rand(length(feedback_scaling),1) < config.mut_rate));
    feedback_scaling(pos) = mutateWeight(feedback_scaling(pos),config);
    offspring.feedback_scaling = reshape(feedback_scaling,size(offspring.feedback_scaling));
    
    feedback_weights = offspring.feedback_weights(:);
    pos = randperm(length(feedback_weights),sum(rand(length(feedback_weights),1) < config.mut_rate));
       
    feedback_weights(pos) = mutateWeight(feedback_weights(pos),config);
    offspring.feedback_weights = reshape(feedback_weights,size(offspring.feedback_weights));
end
end

function value = mutateWeight(value,config)

switch(config.mutate_type)
    case 'gaussian'
        value = value-randn(size(value))*0.15;
        
    case 'uniform'
        value = 2*rand(size(value))-1;
end
end