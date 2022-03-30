%% Infection phase
function loser = recombMinRoRMTSplus(winner,loser,config)

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

% params - W_scaling
W= winner.init_W_scaling(:);
L = loser.init_W_scaling(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.init_W_scaling = reshape(L,size(loser.init_W_scaling));

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
W = winner.W(:);
L = loser.W(:);
indices = find(~loser.test_mask{1});
pos = randperm(length(indices),sum(rand(length(indices),1) < config.rec_rate)); 
L(indices(pos)) = W(indices(pos));
loser.W = reshape(L,size(loser.W));

% delay Weights
W = winner.W_delay(:);
L = loser.W_delay(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate)); 
L(pos) = W(pos);
loser.W_delay = reshape(L,size(loser.W_delay));

% subres connecting weights
if isempty(config.weight_fcn)
    W = winner.W(:);
    L = loser.W(:);
    indices = find(loser.test_mask{2});
    pos = randperm(length(indices),sum(rand(length(indices),1) < config.rec_rate));
    L(indices(pos)) = W(indices(pos));
    loser.W = reshape(L,size(loser.W));
end

% mutate activ fcns
if config.multi_activ
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

%update params
loser = updateParams(loser,config);

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

% MTS
if config.per_node_time_scale % each sub reservoir has the same time update
    individual.update_cycle = individual.init_update_cycle;
else    
    individual.update_cycle(:,1:end_pos) = individual.init_update_cycle(1);
    individual.update_cycle(1:end_pos,:) = individual.init_update_cycle(1);
end

% update subres scaling
for i = 2:length(individual.nodes)
    % set subres positions
    start_pos = end_pos+1;
    end_pos = start_pos + individual.nodes(i)-1;
    
    if ~config.multi_leak_rate
        %individual.input_scaling(:,start_pos:end_pos) = individual.init_input_scaling(:,i);
        individual.leak_rate(start_pos:end_pos) = individual.init_leak_rate(i); 
    end
    
    individual.input_scaling(:,start_pos:end_pos) = individual.init_input_scaling(:,i);
    individual.W_scaling(start_pos:end_pos,start_pos:end_pos) = individual.init_W_scaling(i);
     
    % MTS
    if ~config.per_node_time_scale % each sub reservoir has the same time update
        individual.update_cycle(:,start_pos:end_pos) = individual.init_update_cycle(i); 
        individual.update_cycle(start_pos:end_pos,:) = individual.init_update_cycle(i); 
    end
    
end

% indentify any W that are new connections and put W_scaling as identity
individual.W_scaling((individual.W.*individual.test_mask{1}) ~= 0) = 1;
% check zero connections do not have a scaling
individual.W_scaling(logical(individual.test_mask{1}.*((individual.W.*individual.test_mask{1}) == 0))) = 0;

end