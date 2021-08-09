%% Mutation operator used for different reservoir systems
% Details:
% - number of weights mutated is based on mut_rate;
function offspring = mutateMinRoR(offspring,config)

% input scaling update
input_scaling = offspring.init_input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),[-1 1],config);
offspring.init_input_scaling = reshape(input_scaling,size(offspring.init_input_scaling));

% leak rate update
leak_rate = offspring.init_leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0 1],config);
offspring.init_leak_rate = reshape(leak_rate,size(offspring.init_leak_rate));

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
minimum_weight = 0.01; % anything less than this will be set to zero
prob_2_del = config.prune_rate;

% mutate all input weights
if ~config.CPPN_on
    input_weights = offspring.input_weights(:);
    
    pos = randperm(length(input_weights),sum(rand(length(input_weights),1) < config.mut_rate)); % find new positions to change
    input_weights(pos) = mutateWeight(input_weights(pos),[-1 1],config); % mutate any random weights, chance to mutate existing and non-existing weight
    
    indices = find(input_weights); % find outputs in use
    del_pos = randperm(length(indices),sum(rand(length(indices),1) < prob_2_del)); %get weights to delete
    input_weights(indices(del_pos)) = 0; % delete weight
    
    input_weights(input_weights<minimum_weight & input_weights~=0 & input_weights>-minimum_weight) = 0; % remove small weights
    offspring.input_weights = reshape(input_weights,size(offspring.input_weights));
    
    % mutate subres internal weights    
    if isempty(config.weight_fcn) % only evolve if not given some rigid structure
        W = offspring.W(:); 
        % find all intenal units
        indices = find(~offspring.test_mask{1});
        
        pos = randperm(length(indices),sum(rand(length(indices),1) < config.mut_rate));
        W(indices(pos)) = mutateWeight(W(indices(pos)),[-1 1],config);
        
        % select weights to change
        del_pos = randperm(length(indices),sum(rand(length(indices),1) < prob_2_del));
        W(indices(del_pos)) = 0; % delete weight
        
        W(W<minimum_weight & W~=0 & W>-minimum_weight) = 0;
        offspring.W = reshape(W,size(offspring.W));
    end
    
    % mutate subres connecting weights
    W = offspring.W(:);
    % find all intenal units
    indices = find(offspring.test_mask{2});
    pos = randperm(length(indices),sum(rand(length(indices),1) < config.mut_rate_connecting));
    W(indices(pos)) = mutateWeight(W(indices(pos)),[-1 1],config);
    
    % select weights to change
    del_pos = randperm(length(indices),sum(rand(length(indices),1) < prob_2_del));
    W(indices(del_pos)) = 0; % delete weight
    
    W(W<minimum_weight & W~=0 & W>-minimum_weight) = 0;
    offspring.W = reshape(W,size(offspring.W));
end

% random prune
%offspring = pruneParams(offspring,config);

%update params
offspring = updateParams(offspring,config);

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    
    pos = randperm(length(output_weights),sum(rand(length(output_weights),1) < config.mut_rate)); % find new positions to change
    output_weights(pos) = mutateWeight(output_weights(pos),[-config.output_weight_scaler config.output_weight_scaler],config); % mutate any random weights, chance to mutate existing and non-existing weight
    
    indices = find(output_weights); % find outputs in use
    del_pos = randperm(length(indices),sum(rand(length(indices),1) < prob_2_del)); %get weights to delete
    output_weights(indices(del_pos)) = 0; % delete weight
   
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
    
    % feedback weights
    feedback_weights = offspring.feedback_weights(:);
    pos = randperm(length(feedback_weights),sum(rand(length(feedback_weights),1) < config.mut_rate)); %get weights to delete
    feedback_weights(pos) = mutateWeight(feedback_weights(pos),[-1 1],config); % mutate any random weights, chance to mutate existing and non-existing weight
    
    indices = find(feedback_weights); % find outputs in use
    del_pos = randperm(length(indices),sum(rand(length(indices),1) < prob_2_del)); %get weights to delete
    feedback_weights(indices(del_pos)) = 0; % delete weight

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

%if rand < 0.25
    % prune input
    input_weights = individual.input_weights(:);
    indices = find(input_weights);
    pos =  randperm(length(indices),sum(rand(length(indices),1) < config.prune_rate));
    input_weights(indices(pos)) = 0;
    individual.input_weights = reshape(input_weights,size(individual.input_weights));
%else
    % prune W
    W = individual.W(:);
    indices = find(W);
    pos =  randperm(length(indices),sum(rand(length(indices),1) < config.prune_rate));
    W(indices(pos)) = 0;
    individual.W = reshape(W,size(individual.W));
%end

end

function individual = updateParams(individual,config)

% end of first subres
end_pos = individual.nodes(1);

% intialise W and W_scale
individual.W_scaling(1:end_pos,1:end_pos) = individual.init_W_scaling(1);
individual.input_scaling(:,1:end_pos) = individual.init_input_scaling(:,1);

if config.multi_leak_rate % assign leak rates for each node, if used
    individual.leak_rate = individual.init_leak_rate;
else
    individual.leak_rate(1:end_pos) = individual.init_leak_rate(1);
end

% update subres scaling
for i = 2:length(individual.nodes)
    % set subres positions
    start_pos = end_pos+1;
    end_pos = start_pos + individual.nodes(i)-1;
    
    if ~config.multi_leak_rate
        individual.input_scaling(:,start_pos:end_pos) = individual.init_input_scaling(:,i);
    end
    individual.W_scaling(start_pos:end_pos,start_pos:end_pos) = individual.init_W_scaling(i);
    individual.leak_rate(start_pos:end_pos) = individual.init_leak_rate(i);
end

% indentify any W that are new connections and put W_scaling as identity
individual.W_scaling((individual.W.*individual.test_mask{1}) ~= 0) = 1;
% check zero connections do not have a scaling
individual.W_scaling(logical(individual.test_mask{1}.*((individual.W.*individual.test_mask{1}) == 0))) = 0;

end