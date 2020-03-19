%% Mutation operator used for magnetic film

function offspring = mutateMM(offspring,config)

input_scaling = offspring.input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = input_scaling(pos) + randn*0.15;%2*rand(length(pos),1)-1;
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = leak_rate(pos) + randn*0.15; %rand(length(pos),1);
offspring.leak_rate = reshape(leak_rate,size(offspring.leak_rate));

% vampire params
damping = offspring.damping(:);
pos = randperm(length(damping),sum(rand(length(damping),1) < config.mut_rate));
damping(pos) = mutateWeight(damping(pos),config);
offspring.damping = reshape(damping,size(offspring.damping));

%% fernandos code
% if config.damping_parameter == 'dynamic'
%     if rand < config.mut_rate
%         offspring.damping = offspring.damping - randn*0.15;
%         
%         if offspring.damping < 0
%             offspring.damping = abs(offspring.damping);
%         end
%     end
% end

if config.anisotropy_parameter == 'dynamic'
    if rand < config.mut_rate
        %if rand < 0.5
        %    offspring.anisotropy = 1e-25 + (1e-24-1e-25)*rand;
        %else
        %offspring.anisotropy = 1e-24 + (1e-23-1e-24)*rand;
        offspring.anisotropy = offspring.anisotropy + normrnd(0,1e-23);
        %end
    end
end

if config.temperature_parameter == 'dynamic'
    if rand < config.mut_rate
        offspring.temperature = offspring.temperature + normrnd(0,25);
        
        if offspring.temperature < 0
            offspring.temperature = 0;
        end
    end
end

if config.exchange_parameter == 'dynamic'
    if rand < config.mut_rate
        %offspring.exchange = 1e-21 + (10e-21-1e-21)*rand;
        offspring.exchange = offspring.exchange + normrnd(0,1e-21);
    end
end

if config.magmoment_parameter == 'dynamic'
    if rand < config.mut_rate
        %offspring.magmoment = 0.5 + (5-0.5)*rand;
        offspring.magmoment = offspring.magmoment - randn;
    end
end

if config.applied_field_strength  == 'dynamic'
    if rand < config.mut_rate
        offspring.applied_field_strength = offspring.applied_field_strength - randn*0.15;
    end
end


%% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
    input_weights = offspring.input_weights{i}(:);
    pos =  randperm(length(input_weights),ceil(config.mut_rate*length(input_weights)));
    for n = 1:length(pos)
        input_weights(pos(n)) = input_weights(pos(n)) - randn*0.15;
    end
    offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
    
    % width of inputs
    input_widths = offspring.input_widths{i}(:);
    pos =  randperm(length(input_widths),sum(rand(length(input_widths),1) < config.mut_rate));
    input_widths(pos) = ceil(abs(randn(length(pos),1)));
    input_widths(input_widths > round(sqrt(offspring.nodes(i))/4)) = round(sqrt(offspring.nodes(i))/4);% cap at 1/4 size of space
    offspring.input_widths{i} = reshape(input_widths,size(offspring.input_widths{i}));
    
    % hidden weights
    for j = 1:config.num_reservoirs
        switch(offspring.architecture)
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
