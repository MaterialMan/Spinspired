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

function offspring = mutateHeterotic(offspring,config)

% W scaling
W_scaling = offspring.W_scaling(:);
pos = randperm(length(W_scaling),sum(rand(length(W_scaling),1) < config.mut_rate));
W_scaling(pos) = mutateWeight(W_scaling(pos),[-1 1],config);
offspring.W_scaling = reshape(W_scaling,size(offspring.W_scaling));

leak_rate = offspring.leak_rate(:);
pos =  randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0 1],config);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

% W switch
if isfield(offspring,'W_switch')
    if ~config.RoR_structure
    W_switch = offspring.W_switch(:);
    pos = randperm(length(W_switch),sum(rand(length(W_switch),1) < config.mut_rate));
    W_switch(pos) = round(mutateWeight(W_switch(pos),[0 1],config));
    offspring.W_switch = reshape(W_switch,size(offspring.W_switch));
    end
end

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    % mutate indv res
    offspring.res{i} = offspring.subres_config{i}.mutFcn(offspring.res{i},offspring.subres_config{i});
           
    % mutate connecting weights
    for j = 1:config.num_reservoirs
        switch(config.architecture)
            case 'pipeline'
                if j == i+1
                    W = offspring.W{i,j}(:);
                    % select weights to change
                    pos = randperm(length(W),sum(rand(length(W),1) < config.mut_rate));
                    W(pos) = mutateWeight(W(pos),[-1 1],config);
                    offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
                end
                
            case 'RoR'
                W = offspring.W{i,j}(:);
                % select weights to change
                pos = randperm(length(W),sum(rand(length(W),1) < config.mut_rate));
                W(pos) = mutateWeight(W(pos),[-1 1],config);
                offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
        end
        offspring.connectivity(i,j) = nnz(offspring.W{i,j})/offspring.total_units.^2;        
    end
end


% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos = randperm(length(output_weights),sum(rand(length(output_weights),1) < config.mut_rate));   
    output_weights(pos) = mutateWeight(output_weights(pos),[-config.output_weight_scaler config.output_weight_scaler],config);
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end

% mutate feedback weights
if config.evolve_feedback_weights
    % feedback scaling
    feedback_scaling = offspring.feedback_scaling(:);
    pos =  randperm(length(feedback_scaling),sum(rand(length(feedback_scaling),1) < config.mut_rate));
    feedback_scaling(pos) = mutateWeight(feedback_scaling(pos),[-1 1],config);
    offspring.feedback_scaling = reshape(feedback_scaling,size(offspring.feedback_scaling));
    
    feedback_weights = offspring.feedback_weights(:);
    pos = randperm(length(feedback_weights),sum(rand(length(feedback_weights),1) < config.mut_rate));
    feedback_weights(pos) = mutateWeight(feedback_weights(pos),[-1 1],config);
    offspring.feedback_weights = reshape(feedback_weights,size(offspring.feedback_weights));
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
