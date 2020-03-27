%% Infection phase - need to update
function loser = recombMM(winner,loser,config)

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

%% magnet params
W= winner.damping(:);
L = loser.damping(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.damping = reshape(L,size(loser.damping));

W= winner.anisotropy(:);
L = loser.anisotropy(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.anisotropy = reshape(L,size(loser.anisotropy));

W= winner.temperature(:);
L = loser.temperature(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.temperature = reshape(L,size(loser.temperature));

W= winner.exchange(:);
L = loser.exchange(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.exchange = reshape(L,size(loser.exchange));

W= winner.magmoment(:);
L = loser.magmoment(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.magmoment = reshape(L,size(loser.magmoment));

W= winner.applied_field_strength (:);
L = loser.applied_field_strength (:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.applied_field_strength = reshape(L,size(loser.applied_field_strength ));

%% inputs
for i = 1:config.num_reservoirs
    
    % input weights
    W= winner.input_weights{i}(:);
    L = loser.input_weights{i}(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));    %sum(rand(length(L),1) < config.rec_rate)     
    L(pos) = W(pos);
    loser.input_weights{i} = reshape(L,size(loser.input_weights{i}));
 
        % input widths
    W= winner.input_widths{i}(:);
    L = loser.input_widths{i}(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.input_widths{i} = reshape(L,size(loser.input_widths{i}));
    
    % inner weights
    for j = 1:config.num_reservoirs
        W= winner.W{i,j}(:);
        L = loser.W{i,j}(:);
        pos = randperm(length(L),ceil(config.rec_rate*length(L)));
        L(pos) = W(pos);
        loser.W{i,j} = reshape(L,size(loser.W{i,j}));
    end
    
end

% for output weights
if config.evolve_output_weights
    W= winner.output_weights(:);
    L = loser.output_weights(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.output_weights = reshape(L,size(loser.output_weights));
end