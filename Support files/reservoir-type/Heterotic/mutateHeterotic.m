%% mutate_ReservoirName_.m
% Template function to mutate the offspring reservoir. Use this as a guide when
% creating a new reservoir.
%
% How this function looks at the end depends on the reservoir. However,
% everything below is typically needed to work with all master scripts.
%
% This is called by the @config.mutateFcn pointer.
%
% Additional details:
% - Number of weights mutated is based on mut_rate; 50% chance to change existing weight or remove it

function offspring = mutateHeterotic(offspring,config)

% W scaling
W_scaling = offspring.W_scaling(:);
pos = randperm(length(W_scaling),sum(rand(length(W_scaling),1) < config.mut_rate));
W_scaling(pos) = mutateWeight(W_scaling(pos),config);%2*rand(length(pos),1);
offspring.W_scaling = reshape(W_scaling,size(offspring.W_scaling));

leak_rate = offspring.leak_rate(:);
pos =  randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),config);%rand(length(pos),1);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    % mutate indv res
    offspring.res{i} = offspring.subres_config{i}.mutFcn(offspring.res{i},offspring.subres_config{i});
           
    % hidden weights
    for j = 1:config.num_reservoirs
        switch(config.architecture)
            case 'pipeline'
                if j == i+1
                    W = offspring.W{i,j}(:);
                    % select weights to change
                    pos =  randperm(length(W),ceil(config.mut_rate*length(W)));
                    for n = 1:length(pos)
                        W(pos(n)) = mutateWeight(W(pos(n)),config);
                    end
                    offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
                end
        end
        
        offspring.connectivity(i,j) = nnz(offspring.W{i,j})/offspring.total_units.^2;
        
    end
end


% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos =  randperm(length(output_weights),ceil(config.mut_rate*length(output_weights)));
    
    for n = 1:length(pos)
      output_weights(pos(n)) =  mutateWeight(output_weights(pos(n)),config);
    end
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end
end

function value = mutateWeight(value,config)

switch(config.mutate_type)
    case 'gaussian'
        value = value-randn*0.15;
        
    case 'uniform'
        value = 2*rand-1;
end
end
