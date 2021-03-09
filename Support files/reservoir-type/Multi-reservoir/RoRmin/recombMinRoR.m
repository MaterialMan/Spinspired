%% Infection phase
function loser = recombMinRoR(winner,loser,config)

% params - input_scaling, leak_rate,
W= winner.input_scaling(:);
L = loser.input_scaling(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.input_scaling = reshape(L,size(loser.input_scaling));

W= winner.leak_rate(:);
L = loser.leak_rate(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.leak_rate = reshape(L,size(loser.leak_rate));

% params - W_scaling
W= winner.W_scaling(:);
L = loser.W_scaling(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.W_scaling = reshape(L,size(loser.W_scaling));

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
indices = find(~loser.test_mask);
pos = randperm(length(indices),sum(rand(length(indices),1) < config.rec_rate)); 
L(indices(pos)) = W(indices(pos));
loser.W = reshape(L,size(loser.W));

% subres connecting weights
W = winner.W(:);
L = loser.W(:);
indices = find(loser.test_mask);
pos = randperm(length(indices),sum(rand(length(indices),1) < config.rec_rate)); 
L(indices(pos)) = W(indices(pos));
loser.W = reshape(L,size(loser.W));

% mutate activ fcns
if config.multi_activ
    W= winner.activ_Fcn_indx;
    L = loser.activ_Fcn_indx;
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.activ_Fcn_indx = reshape(L,size(loser.activ_Fcn));
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