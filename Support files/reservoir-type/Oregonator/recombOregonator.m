%% recomb_ReservoirName_.m
% Template function to recombine/infect an individual to create the new offspring reservoir. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
%
% This is called by the @config.recombFcn pointer.

function loser = recombOregonator(winner,loser,config)

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

% adjust input timing
W= winner.time_period(:);
L = loser.time_period(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.time_period = reshape(L,size(loser.time_period));

W= winner.input_widths(:);
L = loser.input_widths(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.input_widths = reshape(L,size(loser.input_widths));

W= winner.input_length(:);
L = loser.input_length(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));         
L(pos) = W(pos);
loser.input_length = reshape(L,size(loser.input_length));

% cycle through sub-reservoirs
for i = 1:config.num_reservoirs
    
     % check input period
    if loser.input_length(i) > loser.time_period(i)
        loser.input_length(i) = loser.time_period(i);
    end
    
    % input weights
    W= winner.input_weights{i}(:);
    L = loser.input_weights{i}(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));     %sum(rand(length(L),1) < config.rec_rate)     
    L(pos) = W(pos);
    loser.input_weights{i} = reshape(L,size(loser.input_weights{i}));
       
    % mix bit and phi matrices
    W_bit= winner.bitmatrix{i}(:);
    L_bit = loser.bitmatrix{i}(:);
    W_phi= winner.phimatrix{i}(:);
    L_phi = loser.phimatrix{i}(:);
    
    pos = randperm(length(L_bit),sum(rand(length(L_bit),1) < config.rec_rate));     %sum(rand(length(L),1) < config.rec_rate)     
    L_bit(pos) = W_bit(pos);
    L_phi(pos) = W_phi(pos);
    
    loser.bitmatrix{i} = reshape(L_bit,size(loser.bitmatrix{i}));
    loser.phimatrix{i} = reshape(L_phi,size(loser.phimatrix{i}));
            
    
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
