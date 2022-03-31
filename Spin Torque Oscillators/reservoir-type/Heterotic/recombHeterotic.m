%% recomb_ReservoirName_.m
% Template function to recombine/infect an individual to create the new offspring reservoir. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
%
% This is called by the @config.recombFcn pointer.

function loser = recombHeterotic(winner,loser,config)

swapped_sub_res = 0;

% params - W_scaling
W= winner.W_scaling(:);
L = loser.W_scaling(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.W_scaling = reshape(L,size(loser.W_scaling));

W= winner.leak_rate(:);
L = loser.leak_rate(:);
pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
L(pos) = W(pos);
loser.leak_rate = reshape(L,size(loser.leak_rate));

% swap entire subres
%if rand > 0.5
% W= winner.res;
% L = loser.res;
% pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
% L(pos) = W(pos);
% loser.res = reshape(L,size(loser.res));
% swapped_sub_res = 1;
%end
 
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

% cycle through sub-reservoirs
for i = 1:config.num_reservoirs
    
    % recomb indv res
    if ~swapped_sub_res
        loser.res{i} = loser.subres_config{i}.recFcn(winner.res{i},loser.res{i},loser.subres_config{i});
    end
    
    % inner weights
    for j = 1:config.num_reservoirs
        switch(config.architecture)
            case 'pipeline'
                if j == i+1
                    W= winner.W{i,j}(:);
                    L = loser.W{i,j}(:);
                    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
                    L(pos) = W(pos);
                    loser.W{i,j} = reshape(L,size(loser.W{i,j}));
                end
            case 'RoR'
                W= winner.W{i,j}(:);
                L = loser.W{i,j}(:);
                pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));
                L(pos) = W(pos);
                loser.W{i,j} = reshape(L,size(loser.W{i,j}));
        end
    end
    
end

% for output weights
if config.evolve_output_weights
    W= winner.output_weights(:);
    L = loser.output_weights(:);
    pos = randperm(length(L),sum(rand(length(L),1) < config.rec_rate));        
    L(pos) = W(pos);
    loser.output_weights = reshape(L,size(loser.output_weights));
end
