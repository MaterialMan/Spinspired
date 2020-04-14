function loser = recombDL(winner,loser,config)

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

% cycle through sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
    W= winner.input_weights{i}(:);
    L = loser.input_weights{i}(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.input_weights{i} = reshape(L,size(loser.input_weights{i}));
    
    % reservoir parameters
    W = [winner.eta(i) winner.gamma(i) winner.p(i)];
    L = [loser.eta(i) loser.gamma(i) loser.p(i)];
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    
    loser.eta(i) = L(1);
    loser.gamma(i) = L(2);
    loser.p(i) = L(3);
       
end

% for output weights
if config.evolve_output_weights
    W= winner.output_weights(:);
    L = loser.output_weights(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
    L(pos) = W(pos);
    loser.output_weights = reshape(L,size(loser.output_weights));
end







