%% Mutation operator used for different reservoir systems
% Details:
% - number of weights mutated is based on mut_rate; 
function offspring = mutateRestrictedRoR(offspring,config)

% input scaling update
input_scaling = offspring.init_input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),[-1 1],config);
offspring.init_input_scaling = reshape(input_scaling,size(offspring.init_input_scaling));

% leak rate update
leak_rate = offspring.init_leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0 1],config);
offspring.init_input_scaling = reshape(leak_rate,size(offspring.init_leak_rate));

% W scaling update
W_scaling = offspring.init_W_scaling(:);
pos = randperm(length(W_scaling),sum(rand(length(W_scaling),1) < config.mut_rate));
W_scaling(pos) = mutateWeight(W_scaling(pos),[-1 1],config);
offspring.init_W_scaling = reshape(W_scaling,size(offspring.init_W_scaling));

% mutate activ fcns
if config.multi_activ
    activFcn = offspring.activ_Fcn_indx;
    pos =  randperm(length(activFcn),sum(rand(length(activFcn),1) < config.mut_rate));
    activFcn(pos) = randi([1 length(config.activ_list)],length(pos),1);
    offspring.activ_Fcn_indx = activFcn;
end
        
% W switch
if isfield(offspring,'W_switch')
    if ~config.RoR_structure
        W_switch = offspring.W_switch(:);
        pos = randperm(length(W_switch),sum(rand(length(W_switch),1) < config.mut_rate));
        W_switch(pos) = round(mutateWeight(W_switch(pos),[0 1],config));
        offspring.W_switch = reshape(W_switch,size(offspring.W_switch));
    end
end

%define minimum weight
minimum_weight = 0.1; % anything less than this will be set to zero

% mutate all input weights
input_weights = offspring.input_weights(:);
pos = randperm(length(input_weights),sum(rand(length(input_weights),1) < config.mut_rate));
input_weights(pos) = mutateWeight(input_weights(pos),[-1 1],config);
input_weights(input_weights<minimum_weight & input_weights~=0 & input_weights>-minimum_weight) = 0;
offspring.input_weights = reshape(input_weights,size(offspring.input_weights));

% restricted_input_mask
restricted_input_mask = offspring.restricted_input_mask(:);
f = find(restricted_input_mask); % find inputs in use
pos = randperm(length(f),sum(rand(length(f),1) < config.mut_rate));
restricted_input_mask(f(pos)) = 0;
restricted_input_mask(randperm(length(restricted_input_mask),length(pos))) = 1;
offspring.restricted_input_mask = reshape(restricted_input_mask,size(offspring.restricted_input_mask));

% restricted_output_mask
restricted_output_mask = offspring.restricted_output_mask(:);
f = find(restricted_output_mask); % find inputs in use
pos = randperm(length(f),sum(rand(length(f),1) < config.mut_rate));
restricted_output_mask(f(pos)) = 0;
restricted_output_mask(randperm(length(restricted_output_mask),length(pos))) = 1;
offspring.restricted_output_mask = reshape(restricted_output_mask,size(offspring.restricted_output_mask));

% mutate subres internal weights
W = offspring.W(:);
% find all intenal units
indices = find(~offspring.test_mask);
% select weights to change
pos = randperm(length(indices),sum(rand(length(indices),1) < config.mut_rate));
W(indices(pos)) = mutateWeight(W(indices(pos)),[-1 1],config);
W(W<minimum_weight & W~=0 & W>-minimum_weight) = 0;
offspring.W = reshape(W,size(offspring.W));

% mutate subres connecting weights
W = offspring.W(:);
% find all intenal units
indices = find(offspring.test_mask);
% select weights to change
pos = randperm(length(indices),sum(rand(length(indices),1) < config.mut_rate_connecting));
W(indices(pos)) = mutateWeight(W(indices(pos)),[-1 1],config);
W(W<minimum_weight & W~=0 & W>-minimum_weight) = 0;
offspring.W = reshape(W,size(offspring.W));

% random prune
offspring = pruneParams(offspring,config);

%update params
offspring = updateParams(offspring);

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos = randperm(length(output_weights),sum(rand(length(output_weights),1) < config.mut_rate));
    output_weights(pos) = mutateWeight(output_weights(pos),[-config.output_weight_scaler config.output_weight_scaler],config);
    output_weights(output_weights<minimum_weight & output_weights~=0 & output_weights>-minimum_weight) = 0;
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
    feedback_weights(feedback_weights<minimum_weight & feedback_weights~=0 & feedback_weights>-minimum_weight) = 0;
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


function individual = pruneParams(individual,config)

if rand < 0.25
    % prune input
    input_weights = individual.input_weights(:);
    indices = find(input_weights);
    pos =  randperm(length(indices),sum(rand(length(indices),1) < config.prune_rate));
    input_weights(indices(pos)) = 0;
    individual.input_weights = reshape(input_weights,size(individual.input_weights));
else
    % prune W
    W = individual.W(:);
    indices = find(W);
    pos =  randperm(length(indices),sum(rand(length(indices),1) < config.prune_rate));
    W(indices(pos)) = 0;
    individual.W = reshape(W,size(individual.W));
end

end

function individual = updateParams(individual)

% end of first subres
end_pos = individual.nodes(1);

% intialise W and W_scale
individual.W_scaling(1:end_pos,1:end_pos) = individual.init_W_scaling(1);
individual.input_scaling(:,1:end_pos) = individual.init_input_scaling(:,1);
individual.leak_rate(1:end_pos) = individual.init_leak_rate(1);

% update subres scaling
for i = 2:length(individual.nodes)
    
    start_pos = end_pos+1;
    end_pos = start_pos + individual.nodes(i)-1;
    
    individual.input_scaling(:,start_pos:end_pos) = individual.init_input_scaling(:,i);
    individual.W_scaling(start_pos:end_pos,start_pos:end_pos) = individual.init_W_scaling(i);
    individual.leak_rate(start_pos:end_pos) = individual.init_leak_rate(i);        
end

% indentify any W that are new connections and put W_scaling as identity
individual.W_scaling((individual.W.*individual.test_mask) ~= 0) = 1;
% check zero connections do not have a scaling
individual.W_scaling(logical(individual.test_mask.*((individual.W.*individual.test_mask) == 0))) = 0;

end