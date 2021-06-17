%% Mutation operator used for different reservoir systems
% Details:
% - number of weights mutated is based on mut_rate; 
function offspring = mutateRoR(offspring,config)

% params - input scaling and leak rate
input_scaling = offspring.input_scaling(:);
pos = randperm(length(input_scaling),sum(rand(length(input_scaling),1) < config.mut_rate));
input_scaling(pos) = mutateWeight(input_scaling(pos),[-1 1],config);
offspring.input_scaling = reshape(input_scaling,size(offspring.input_scaling));

leak_rate = offspring.leak_rate(:);
pos = randperm(length(leak_rate),sum(rand(length(leak_rate),1) < config.mut_rate));
leak_rate(pos) = mutateWeight(leak_rate(pos),[0 1],config);
offspring.input_scaling = reshape(leak_rate,size(offspring.leak_rate));

% W scaling
W_scaling = offspring.W_scaling(:);
pos = randperm(length(W_scaling),sum(rand(length(W_scaling),1) < config.mut_rate));
W_scaling(pos) = mutateWeight(W_scaling(pos),[-1 1],config);
offspring.W_scaling = reshape(W_scaling,size(offspring.W_scaling));

% W switch
if isfield(offspring,'W_switch')
    if ~config.RoR_structure
        W_switch = offspring.W_switch(:);
        pos = randperm(length(W_switch),sum(rand(length(W_switch),1) < config.mut_rate));
        W_switch(pos) = round(mutateWeight(W_switch(pos),[0 1],config));
        offspring.W_switch = reshape(W_switch,size(offspring.W_switch));
    end
end

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    %anything below this becomes zero
    minimum_weight = 0.1;
    
    % input weights
    input_weights = offspring.input_weights{i}(:);
    pos = randperm(length(input_weights),sum(rand(length(input_weights),1) < config.mut_rate));
    input_weights(pos) = mutateWeight(input_weights(pos),[-1 1],config); 
    input_weights(input_weights<minimum_weight & input_weights~=0 & input_weights>-minimum_weight) = 0;
    offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
        
    % hidden weights
    for j = 1:config.num_reservoirs
        % only mutate one half of matrix if undirected weights in use
        if (config.undirected_ensemble && i ~= j) || (config.undirected && i == j)
            W= triu(offspring.W{i,j});
            f = find(W);
             pos = randperm(length(f),sum(rand(length(f),1) < config.mut_rate));
            W(f(pos)) = mutateWeight(W(f(pos)),[-1 1],config);
            W = triu(W)+triu(W,1)'; % copy top-half to lower-half
            offspring.W{i,j} = W;
        else
            W = offspring.W{i,j}(:);
            % select weights to change
            pos = randperm(length(W),sum(rand(length(W),1) < config.mut_rate));
            W(pos) = mutateWeight(W(pos),[-1 1],config);
            
            W(W<minimum_weight & W~=0 & W>-minimum_weight) = 0;
            offspring.W{i,j} = reshape(W,size(offspring.W{i,j}));
        end
               
        offspring.connectivity(i,j) = nnz(offspring.W{i,j})/offspring.total_units.^2;
    end
    
    % mutate activ fcns
    if config.multi_activ
        activFcn = offspring.activ_Fcn(i,1:offspring.nodes(i));
        pos =  randperm(length(activFcn),sum(rand(length(activFcn),1) < config.mut_rate));
        activFcn(pos) = {config.activ_list{randi([1 length(config.activ_list)],length(pos),1)}};
        offspring.activ_Fcn(i,1:offspring.nodes(i)) = reshape(activFcn,size(offspring.activ_Fcn(i,1:offspring.nodes(i))));
    else
        activFcn = offspring.activ_Fcn;
        pos =  randperm(length(activFcn),sum(rand(length(activFcn),1) < config.mut_rate));
        activFcn(pos) = {config.activ_list{randi([1 length(config.activ_list)],length(pos),1)}};
        offspring.activ_Fcn = reshape(activFcn,size(offspring.activ_Fcn));
    end
    
    if config.iir_filter_on
        iir_feedfoward = offspring.iir_weights{i,1}(:,1);
        pos = randperm(length(iir_feedfoward),sum(rand(length(iir_feedfoward),1) < config.mut_rate));
        w_0 = mutateWeight(iir_feedfoward(pos),[-1 1],config);        
        alpha = sin(w_0).*sinh((log(2)./2) * (3*rand) * (w_0./(sin(w_0))));
        offspring.iir_weights{i,1}(pos,:) = alpha .* [1 0 -1]; 
                
        %offspring.iir_weights{i,1} = reshape(iir_feedfoward,size(offspring.iir_weights{i,1}));
        
        %iir_feedback = offspring.iir_weights{i,2}(:,1);
        %pos =  randperm(length(iir_feedback),ceil(config.mut_rate*length(iir_feedback)));
        %iir_feedback(pos) = mutateWeight(iir_feedback(pos),config);
        
        offspring.iir_weights{i,2}(pos,:) = [1+alpha -2*cos(w_0) 1-alpha];%reshape(iir_feedback,size(offspring.iir_weights{i,2}));        
        
    end
end

% mutate output weights
% if config.evolve_output_weights
%     output_weights = offspring.output_weights(:);
%     pos = randperm(length(output_weights),sum(rand(length(output_weights),1) < config.mut_rate));   
%     output_weights(pos) = mutateWeight(output_weights(pos),[-config.output_weight_scaler config.output_weight_scaler],config);
%     output_weights(output_weights<minimum_weight & output_weights~=0 & output_weights>-minimum_weight) = 0;
%     offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
% end

% mutate feedback weights
if config.evolve_feedback_weights
    % feedback scaling
    feedback_scaling = offspring.feedback_scaling(:);
    pos =  randperm(length(feedback_scaling),sum(rand(length(feedback_scaling),1) < config.mut_rate));
    feedback_scaling(pos) = mutateWeight(feedback_scaling(pos),[-1 1],config);
    offspring.feedback_scaling = reshape(feedback_scaling,size(offspring.feedback_scaling));
    
    feedback_weights = offspring.feedback_weights(:);
    pos = randperm(length(feedback_weights),sum(rand(length(feedback_weights),1) < config.mut_rate));
    feedback_weights(pos) = mutateWeight(feedback_weights(pos),[-1 1],config);
    feedback_weights(feedback_weights<minimum_weight & feedback_weights~=0 & feedback_weights>-minimum_weight) = 0;
    offspring.feedback_weights = reshape(feedback_weights,size(offspring.feedback_weights));
end
end

function value = mutateWeight(value,range,config)

if range(1)~=range(2)
    switch(config.mutate_type)
        case 'gaussian'
            for i = 1:length(value)
                flag = 1;
                while(flag)
                    t_value = value(i) + (range(1) + (range(2)-range(1))*randn);
                    
                    % check within range
                    if (t_value <= range(2)) && (t_value >= range(1))
                        flag = 0;
                    end
                end
                value(i) = t_value;
            end
        case 'uniform'
            value = range(1) + (range(2)-range(1))*rand;
    end
end
end