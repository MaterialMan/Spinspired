%% Infection phase
function loser = recombMinRoRMTS(winner,loser,config)

% params - input_scaling, leak_rate,
W= winner.init_input_scaling(:);
L = loser.init_input_scaling(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.init_input_scaling = reshape(L,size(loser.init_input_scaling));

W= winner.init_leak_rate(:);
L = loser.init_leak_rate(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.init_leak_rate = reshape(L,size(loser.init_leak_rate));

W= winner.init_update_cycle(:);
L = loser.init_update_cycle(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.init_update_cycle = reshape(L,size(loser.init_update_cycle));

W= winner.init_max_update_cycle(:);
L = loser.init_max_update_cycle(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.init_max_update_cycle = reshape(L,size(loser.init_max_update_cycle));

W= winner.init_input_delay(:);
L = loser.init_input_delay(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.init_input_delay = reshape(L,size(loser.init_input_delay));

% params - W_scaling
W= winner.init_W_scaling(:);
L = loser.init_W_scaling(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.init_W_scaling = reshape(L,size(loser.init_W_scaling));

% params - interpolation_length
W= winner.interpolation_length(:);
L = loser.interpolation_length(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.interpolation_length = reshape(L,size(loser.interpolation_length));

% params - quatisation 2^m
if config.quantized_state(1) > 0
W= winner.m(:);
L = loser.m(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.m = reshape(L,size(loser.m));
end

% W switch
if isfield(winner,'W_switch')
    if ~config.RoR_structure
        W= winner.W_switch(:);
        L = loser.W_switch(:);
        pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
        L(pos) = W(pos);
        loser.W_switch = reshape(L,size(loser.W_switch));
    end
end

% input weights
W= winner.input_weights(:);
L = loser.input_weights(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));     %sum(rand(length(L),1) < config.rec_rate)
L(pos) = W(pos);
loser.input_weights = reshape(L,size(loser.input_weights));

% subres internal weights
if ~config.dipole_fields
    W = winner.W(:);
    L = loser.W(:);
    indices = find(~loser.test_mask{1});
    pos = randperm(length(indices),sum(rand(length(indices),1) < config.rec_rate));
    L(indices(pos)) = W(indices(pos));
    loser.W = reshape(L,size(loser.W));   
end

% subres connecting weights
if ~contains(config.internal_weight_initialisation,["graph","weight_fcn","dipole_fields"]) && ~isempty(config.RoR_structure)%isempty(config.weight_fcn)
    W = winner.W(:);
    L = loser.W(:);
    indices = find(loser.test_mask{2});
    pos = randperm(length(indices),sum(rand(length(indices),1) < config.rec_rate));
    L(indices(pos)) = W(indices(pos));
    loser.W = reshape(L,size(loser.W));
end
    
% mutate activ fcns
if length(config.activ_list) > 1
    W= winner.activ_Fcn_indx;
    L = loser.activ_Fcn_indx;
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.activ_Fcn_indx = reshape(L,size(loser.activ_Fcn_indx));
end


% for output weights
if config.evolve_output_weights
    W= winner.output_weights(:);
    L = loser.output_weights(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.output_weights = reshape(L,size(loser.output_weights));
end

% for feedback weights
if config.evolve_feedback_weights
    % params - W_scaling
    W= winner.feedback_scaling(:);
    L = loser.feedback_scaling(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.feedback_scaling = reshape(L,size(loser.feedback_scaling));
    
    W= winner.feedback_weights(:);
    L = loser.feedback_weights(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.feedback_weights = reshape(L,size(loser.feedback_weights));
end

if config.dipole_fields
    
    % reset flag
    flag = zeros(4,1);
    
    % perim_value
    W= winner.perim_value(:);
    L = loser.perim_value(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    if ~isempty(pos)
        flag(1) = 1;
    end
    loser.perim_value = reshape(L,size(loser.perim_value));
    
    % self_feed_back_loop
    W= winner.self_feed_back_loop(:);
    L = loser.self_feed_back_loop(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    if ~isempty(pos)
        flag(2) = 1;
    end
    loser.self_feed_back_loop = reshape(L,size(loser.self_feed_back_loop));
    
    % dipole_field_value
    W= winner.dipole_scale(:);
    L = loser.dipole_scale(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    if ~isempty(pos)
        flag(3) = 1;
    end
    loser.dipole_scale = reshape(L,size(loser.dipole_scale));
    
else
    flag = 0;
end

%update params
loser = updateParams(loser,config,flag);

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
    %individual.max_update_cycle(1:end_pos) = individual.init_max_update_cycle(1);
end

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
    
    individual.input_scaling(:,start_pos:end_pos) = individual.init_input_scaling(:,i);
    individual.W_scaling(start_pos:end_pos,start_pos:end_pos) = individual.init_W_scaling(i);
    
    % MTS
    if ~config.per_node_time_scale % each sub reservoir has the same time update
        % check new update cycles are less then new max update
        if individual.init_update_cycle(i) > individual.init_max_update_cycle
            individual.init_update_cycle(i) = individual.init_max_update_cycle;
        end
        individual.update_cycle(start_pos:end_pos) = individual.init_update_cycle(i);
        %individual.max_update_cycle(start:end_pos) = individual.init_max_update_cycle;
    end
    
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


end