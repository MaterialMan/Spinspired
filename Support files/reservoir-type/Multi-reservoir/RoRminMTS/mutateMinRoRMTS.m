%% Mutation operator used for different reservoir systems
% Details:
% - number of weights mutated is based on mut_rate;
function offspring = mutateMinRoRMTS(offspring,config)

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

% cycle update
update_cycle = offspring.init_update_cycle(:);
pos = randperm(length(update_cycle),sum(rand(length(update_cycle),1) < config.mut_rate));
update_cycle(pos) = floor(mutateWeight(update_cycle(pos),[0 config.max_update_cycle],config));
del_pos = randperm(length(update_cycle),sum(rand(length(update_cycle),1) < config.prune_rate));
update_cycle(del_pos) = 1; % delete input delay
offspring.init_update_cycle = reshape(update_cycle,size(offspring.init_update_cycle));

% the maximum cycle update
max_update_cycle = offspring.init_max_update_cycle(:);
pos = randperm(length(max_update_cycle),sum(rand(length(max_update_cycle),1) < config.mut_rate));
max_update_cycle(pos) = floor(mutateWeight(max_update_cycle(pos),[0 config.max_update_mutate_max],config));
offspring.init_max_update_cycle = reshape(max_update_cycle,size(offspring.init_max_update_cycle));

% input delay
input_delay = offspring.init_input_delay(:);
pos =  randperm(length(input_delay),sum(rand(length(input_delay),1) < config.mut_rate));
input_delay(pos) = floor(mutateWeight(input_delay(pos),[0,config.max_input_delay],config));
del_pos = randperm(length(input_delay),sum(rand(length(input_delay),1) < config.prune_rate));
input_delay(del_pos) = 0; % delete input delay
offspring.init_input_delay = reshape(input_delay,size(offspring.init_input_delay));

% W scaling update
W_scaling = offspring.init_W_scaling(:);
pos = randperm(length(W_scaling),sum(rand(length(W_scaling),1) < config.mut_rate));
W_scaling(pos) = mutateWeight(W_scaling(pos),[-1 1],config);
offspring.init_W_scaling = reshape(W_scaling,size(offspring.init_W_scaling));

% interpolation_length
interpolation_length = offspring.interpolation_length(:);
pos = randperm(length(interpolation_length),sum(rand(length(interpolation_length),1) < config.mut_rate));
interpolation_length(pos) = round(mutateWeight(interpolation_length(pos),[1 config.max_interpolation_length],config));
offspring.interpolation_length = reshape(interpolation_length,size(offspring.interpolation_length));

% Quantized m (2^m states)
if config.quantized_state(1) > 0
    quatize_m = offspring.m(:);
    pos = randperm(length(quatize_m),sum(rand(length(quatize_m),1) < config.mut_rate));
    quatize_m(pos) = round(mutateWeight(quatize_m(pos),[config.quantized_state(1) config.quantized_state(2)],config));
    offspring.m = reshape(quatize_m,size(offspring.m));
end

% mutate activ fcns
if length(config.activ_list) > 1
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
    if ~config.dipole_fields % only evolve if not given some rigid structure
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
    
    if  ~contains(config.internal_weight_initialisation,["graph","weight_fcn","dipole_fields"]) && ~isempty(config.RoR_structure)%isempty(config.weight_fcn)
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
end

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

% alter parameters that control dipole configuration
if config.dipole_fields
    % reset flag
    flag = zeros(4,1);
    
    % perim_value
    perim_value = offspring.perim_value(:);
    pos = randperm(length(perim_value),sum(rand(length(perim_value),1) < config.mut_rate));
    perim_value(pos) = mutateWeight(perim_value(pos),[-1 1],config);
    %del_pos = randperm(length(perim_value),sum(rand(length(perim_value),1) < prob_2_del));
    %perim_value(del_pos) = 0; % delete weight
    if ~isequal(perim_value,offspring.perim_value(:))
        flag(1) = 1;
    end
    offspring.perim_value = reshape(perim_value,size(offspring.perim_value));
    
    % self_feed_back_loop
    self_feed_back_loop = offspring.self_feed_back_loop(:);
    pos = randperm(length(self_feed_back_loop),sum(rand(length(self_feed_back_loop),1) < config.mut_rate));
    self_feed_back_loop(pos) = mutateWeight(self_feed_back_loop(pos),[-1 1],config);
    %del_pos = randperm(length(self_feed_back_loop),sum(rand(length(self_feed_back_loop),1) < prob_2_del));
    %self_feed_back_loop(del_pos) = 0; % delete weight
    if ~isequal(self_feed_back_loop,offspring.self_feed_back_loop(:))
        flag(2) = 1;
    end
    offspring.self_feed_back_loop = reshape(self_feed_back_loop,size(offspring.self_feed_back_loop));
    
    % dipole_field_value
    dipole_scale= offspring.dipole_scale;
    pos = randperm(length(dipole_scale),sum(rand(length(dipole_scale),1) < config.mut_rate));
    dipole_scale(pos) = mutateWeight(dipole_scale(pos),[-1 1],config);
    if ~isequal(dipole_scale,offspring.dipole_scale)
        flag(3) = 1;
    end
    offspring.dipole_scale = reshape(dipole_scale,size(offspring.dipole_scale));
    
    % num_dipoles
    num_dipoles = offspring.num_dipoles;
    pos = randperm(length(num_dipoles),sum(rand(length(num_dipoles),1) < config.mut_rate));
    num_dipoles(pos) = round(mutateWeight(num_dipoles(pos),[1 config.num_dipoles],config));
    if ~isequal(num_dipoles,offspring.num_dipoles)
        flag(4) = 1;
    end
    offspring.num_dipoles = reshape(num_dipoles,size(offspring.num_dipoles));
    
    
else
    flag = 0;
end

% random prune
%offspring = pruneParams(offspring,config);

%update params
offspring = updateParams(offspring,config,flag);

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

function individual = updateParams(individual,config,flag)

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

% MTS
if config.per_node_time_scale % each sub reservoir has the same time update
    % check new update cycles are less then new max update
    individual.init_update_cycle(individual.init_update_cycle > individual.init_max_update_cycle) = individual.init_max_update_cycle;
    individual.update_cycle = individual.init_update_cycle;
    %individual.max_update_cycle = individual.init_max_update_cycle;
else
    % check new update cycles are less then new max update
    if individual.init_update_cycle(1) > individual.init_max_update_cycle
        individual.init_update_cycle(1) = individual.init_max_update_cycle;
    end
    individual.update_cycle(1:end_pos) = individual.init_update_cycle(1);
    %individual.max_update_cycle(1:end_pos) = individual.init_max_update_cycle(1);
end

% input delay
if config.per_input_delay
    individual.input_delay = individual.init_input_delay;
else
    individual.input_delay(1:end_pos) = individual.init_input_delay(1);
end

% update subres scaling
for i = 2:length(individual.nodes)
    % set subres positions
    start_pos = end_pos+1;
    end_pos = start_pos + individual.nodes(i)-1;
    
    if ~config.multi_leak_rate
        individual.leak_rate(start_pos:end_pos) = individual.init_leak_rate(i);
    end
    individual.W_scaling(start_pos:end_pos,start_pos:end_pos) = individual.init_W_scaling(i);
    individual.input_scaling(:,start_pos:end_pos) = individual.init_input_scaling(:,i);
    
    % MTS
    if ~config.per_node_time_scale % each sub reservoir has the same time update
        % check new update cycles are less then new max update
        if individual.init_update_cycle(i) > individual.init_max_update_cycle
            individual.init_update_cycle(i) = individual.init_max_update_cycle;
        end
        individual.update_cycle(start_pos:end_pos) = individual.init_update_cycle(i);
        %individual.max_update_cycle(start_pos:end_pos) = individual.init_max_update_cycle(i);
    end
    
    % input delay
    if ~config.per_input_delay
        individual.input_delay(start_pos:end_pos) = individual.init_input_delay(i);
    end
end

% indentify any W that are new connections and put W_scaling as identity
individual.W_scaling((individual.W.*individual.test_mask{1}) ~= 0) = 1;
% check zero connections do not have a scaling
individual.W_scaling(logical(individual.test_mask{1}.*((individual.W.*individual.test_mask{1}) == 0))) = 0;

if config.dipole_fields
    
    for res = 1:length(individual.nodes)
        
        if res == 1
            start_pos = 1;
            end_pos = individual.nodes(1);
        else
            % set subres positions
            start_pos = end_pos+1;
            end_pos = start_pos + individual.nodes(res)-1;
        end
        
        % get weights to change
        weights = individual.W(start_pos:end_pos,start_pos:end_pos);
        
        % set strength of dipole fields
        if flag(4) || flag(3)
            dipole_scale = individual.dipole_scale; %0.2
            dipole_field_value = (1./((1:individual.nodes(res)).^3))*dipole_scale;%1./exp((1:config.num_nodes)*1.1)*3;% 1./exp((1:config.num_nodes)*1.1)*3
            
            for n = 1:individual.nodes(res)
                source_node = n;
                for i = 1:individual.num_dipoles
                    [nn,dist] = nearest(config.G{res},source_node,i,'Method','unweighted');
                    weights(sub2ind(size(weights),repmat(source_node,length(nn(dist == i)),1),nn(dist == i))) = dipole_field_value(i);
                    weights(sub2ind(size(weights),nn(dist == i),repmat(source_node,length(nn(dist == i)),1))) = dipole_field_value(i);
                end
            end
            individual.W(start_pos:end_pos,start_pos:end_pos) = weights;
        end
        
        % set feedback weights
        W_scaling = individual.W_scaling(start_pos:end_pos,start_pos:end_pos);
        if flag(2)
            idx = logical(eye(individual.nodes(res),individual.nodes(res)));
            W_scaling(idx) = individual.self_feed_back_loop;
        end
        
        % out nodes
        if flag(1)
            sq = sqrt(individual.nodes(res));
            perim_nodes = [1:sq; 1:sq:individual.nodes(res); sq:sq:individual.nodes(res);
                (1:sq)+(individual.nodes(res)-sq)];
            
            %perim_scaling = -individual.init_W_scaling(res);%individual.perim_value(1);
            %permin_self_loop_scaling = -individual.self_feed_back_loop;%individual.perim_value(2);
            perim_scaling = individual.perim_value(1);
            permin_self_loop_scaling = individual.perim_value(2);
            config.corner_scale = individual.perim_value(3);
            
            for i = 1:individual.nodes(res)
                
                if ismember(i,unique(perim_nodes))
                    % set back connections
                    W_scaling(i+1:end,i) = perim_scaling;
                    W_scaling(i,i+1:end) = perim_scaling;
                    
                    W_scaling(i:-1:1,i) = perim_scaling;
                    W_scaling(i,i:-1:1) = perim_scaling;
                    
                    % assign feeback connection
                    W_scaling(i,i) = permin_self_loop_scaling;
                end
            end
            
            % update corner links - rescale corner connections
            for i = 1:individual.nodes(res)
                
                if i-sqrt(individual.nodes(res))-1 > 0
                    W_scaling(i,i-sqrt(individual.nodes(res))-1) = W_scaling(i,i-sqrt(individual.nodes(res))-1)*config.corner_scale;
                    
                end
                if i-sqrt(individual.nodes(res))+1 > 0
                    W_scaling(i,i-sqrt(individual.nodes(res))+1) = W_scaling(i,i-sqrt(individual.nodes(res))+1)*config.corner_scale;
                    
                end
                if i+sqrt(individual.nodes(res))-1 > 0 && i > 1 && i+sqrt(individual.nodes(res))+1 < individual.nodes(res)
                    W_scaling(i,i+sqrt(individual.nodes(res))-1) = W_scaling(i,i+sqrt(individual.nodes(res))-1)*config.corner_scale;
                    
                end
                if i+sqrt(individual.nodes(res))+1 > 0 && i+sqrt(individual.nodes(res))+1 < individual.nodes(res)
                    W_scaling(i,i+sqrt(individual.nodes(res))+1) = W_scaling(i,i+sqrt(individual.nodes(res))+1)*config.corner_scale;
                end
            end
        end
        individual.W_scaling(start_pos:end_pos,start_pos:end_pos) = W_scaling;
    end
    
end

% if contains(config.dataset, 'mimic')
%     individual.output_weights = eye(individual.nodes);
%     inputs = zeros(2,100);
%     inputs(1,:) = [0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-0.0152088000000000	-0.0220538000000000	-0.0150787000000000	0	0	0	0	0	0	-0.0152312000000000	-0.0409663000000000	0.0707502000000000	-0.0409777000000000	-0.0163733000000000	0	0	0	0	0	-0.0221606000000000	0.0705514000000000	0.304361000000000	0.0701387000000000	-0.0217581000000000	0	0	0	0	0	-0.0151052000000000	-0.0409832000000000	0.0703002000000000	-0.0409814000000000	-0.0162503000000000	0	0	0	0	0	0	-0.0163581000000000	-0.0216611000000000	-0.0162313000000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
%     individual.input_weights = inputs;
% end

end