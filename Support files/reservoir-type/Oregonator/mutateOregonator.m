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

function offspring = mutateOregonator(offspring,config)
     
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

input_widths = offspring.input_widths(:);
pos =  randperm(length(input_widths),sum(rand(length(input_widths),1) < config.mut_rate));
input_widths(pos) = randi([1 config.vesicle_radius-1],length(pos),1);
offspring.input_widths = reshape(input_widths,size(offspring.input_widths));

time_period = offspring.time_period(:);
pos =  randperm(length(time_period),sum(rand(length(time_period),1) < config.mut_rate));
time_period(pos) = randi([1 config.max_time_period],length(pos),1);
offspring.time_period = reshape(time_period,size(offspring.time_period));

input_length = offspring.input_length(:);
pos =  randperm(length(input_length),sum(rand(length(input_length),1) < config.mut_rate));
input_length(pos) = randi([1 offspring.time_period],length(pos),1);
offspring.input_length = reshape(input_length,size(offspring.input_length));

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs

    % check input period
    if offspring.input_length(i) > offspring.time_period(i)
        offspring.input_length(i) = offspring.time_period(i);
    end

    % input weights
    input_weights = offspring.input_weights{i}(:);
    pos = randperm(length(input_weights),sum(rand(length(input_weights),1) < config.mut_rate));
    input_weights(pos) = mutateWeight(input_weights(pos),config); 
    offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
        
    % change bit and phi matrices
    bitmatrix = offspring.bitmatrix{i}(:);
    phimatrix = offspring.phimatrix{i}(:);
    
    pos = randperm(length(bitmatrix),sum(rand(length(bitmatrix),1) < config.mut_rate));
    bitmatrix(pos) = round(rand(length(pos),1)); 
    to_change = bitmatrix(pos) == 1;
    phimatrix(pos(to_change)) = rand(sum(to_change),1)*0.05;%mutateWeight(phimatrix(pos(to_change)),config);

    offspring.bitmatrix{i} = reshape(bitmatrix,size(offspring.bitmatrix{i}));
    offspring.phimatrix{i} = reshape(phimatrix,size(offspring.phimatrix{i}));           
    
end

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