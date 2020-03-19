%% recomb_ReservoirName_.m
% Template function to recombine/infect an individual to create the new offspring reservoir. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
%
% This is called by the @config.recombFcn pointer.

function loser = recombIsing(winner,loser,config)

% params - input_scaling, leak_rate,
W= winner.input_scaling(:);
L = loser.input_scaling(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
L(pos) = W(pos);
loser.input_scaling = reshape(L,size(loser.input_scaling));

W= winner.leak_rate(:);
L = loser.leak_rate(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
L(pos) = W(pos);
loser.leak_rate = reshape(L,size(loser.leak_rate));

% ising param
W= winner.J(:);
L = loser.J(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));
L(pos) = W(pos);
loser.J = reshape(L,size(loser.J));

W= winner.kT(:);
L = loser.kT(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));
L(pos) = W(pos);
loser.kT = reshape(L,size(loser.kT));

W= winner.probSpinUp(:);
L = loser.probSpinUp(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));
L(pos) = W(pos);
loser.probSpinUp = reshape(L,size(loser.probSpinUp));

% cycle through sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
    W= winner.input_weights{i}(:);
    L = loser.input_weights{i}(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.input_weights{i} = reshape(L,size(loser.input_weights{i}));
       
     % input widths
    W= winner.input_widths{i}(:);
    L = loser.input_widths{i}(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.input_widths{i} = reshape(L,size(loser.input_widths{i})); 
    
    W= winner.init_spin{i}(:);
    L = loser.init_spin{i}(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.init_spin{i} = reshape(L,size(loser.init_spin{i}));       
    
end

% for output weights
if config.evolve_output_weights
    W= winner.output_weights(:);
    L = loser.output_weights(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.output_weights = reshape(L,size(loser.output_weights));
end
