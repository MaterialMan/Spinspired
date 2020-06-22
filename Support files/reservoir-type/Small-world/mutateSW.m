%% Mutation operator used for SW reservoir systems
% Details:
% - number of weights mutated is based on mut_rate;
function offspring = mutateSW(offspring,config)

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
W_scaling(pos) = mutateWeight(W_scaling(pos),[0 2],config);
offspring.W_scaling = reshape(W_scaling,size(offspring.W_scaling));

% cycle through all sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
    input_weights = offspring.input_weights{i}(:);
    pos = randperm(length(input_weights),sum(rand(length(input_weights),1) < config.mut_rate));
    input_weights(pos) = mutateWeight(input_weights(pos),[-1 1],config); 
    offspring.input_weights{i} = reshape(input_weights,size(offspring.input_weights{i}));
    
    % hidden weights
    for j = 1:config.num_reservoirs
        
        switch(config.SW_type)
            
            case 'topology'
                
                W = offspring.W{i,j};
                %change base graph
                f = find(adjacency(config.G{1,j}));
                pos = randperm(length(f),sum(rand(length(f),1) < config.mut_rate));
                W(f(pos)) = mutateWeight(W(f(pos)),[-1 1],config);
                offspring.W{i,j} = W;
                
            case 'topology_plus_weights'% must maintain proportion of connections
                W = offspring.W{i,j};  % current graph
                base_W_0 = adjacency(config.G{1,j});
                pos_chng = find(~base_W_0); % non-base weights
                
                w1 = find(W(pos_chng)); %all non-zero non-base weights
                w = w1; %set default non-zero non-base weights
                
                for p = 1:ceil(config.mut_rate*length(w1)) % num to mutate
                    
                    pos(p) = randperm(length(w),1);
                    while(sum(pos(p) == pos(1:p-1)) > 0)
                        pos(p) = randperm(length(w),1);
                    end
                    
                    if round(rand) && config.P_rc < 1
                        % remove random non-zero non-base weight
                        W(pos_chng(w(pos(p)))) = 0;
                        
                        pos2(p) = w(pos(p));
                        while(sum(pos2(p) == w) > 0 || sum(pos2(p) == pos2(1:p-1)) > 0)
                            pos2(p) = randi([1 length(pos_chng)]);
                        end
                        
                        W(pos_chng(pos2(p))) = mutateWeight(W(pos_chng(pos2(p))),[-1 1],config);
                        
                        %check still okay
                        if nnz(offspring.W{i,j}) ~= nnz(W)
                            error('SW not working');
                        end
                    else
                        % change non-zero non-base weight to another value
                        W(pos_chng(w(pos(p)))) = mutateWeight(pos_chng(w(pos(p))),[-1 1],config);
                    end
                    
                end
                
                %check still okay
                if nnz(offspring.W{i,j}) ~= nnz(W)
                    error('SW not working');
                end
                
                offspring.W{i,j} = W;
                
                %change base graph
                f = find(base_W_0);
                pos = randperm(length(f),sum(rand(length(f),1) < config.mut_rate));
                
                % select weights to change
                W(f(pos)) = mutateWeight(W(f(pos)),[-1 1],config);
                offspring.W{i,j} = W;
                
            case 'watts_strogartz'
                
                W = offspring.W{i,j};
                f = find(W);
                pos = randperm(length(f),sum(rand(length(f),1) < config.mut_rate));
                for n = 1:length(pos)
                    % switch
                    if rand < config.P_rc% rewiring probability
                        % find row and col
                        [row,col] = ind2sub(size(W),f(pos(n)));
                        
                        tmp_val1 = W(row,col);
                        tmp_val2 = W(col,row);
                        
                        W(row,col) = 0;
                        W(col,row) = 0;
                        
                        list = 1:length(W);
                        list(list == row) = [];
                        list(list == col) = [];
                        
                        indx = randi([1 length(list)]);
                        
                        W(row,list(indx)) = tmp_val1;
                        W(list(indx),col) = tmp_val2;
                    else
                        W(f(pos(n))) = mutateWeight(W(f(pos(n))),[-1 1],config);
                    end
                end
        end
        
        offspring.connectivity(i,j) = nnz(offspring.W{i,j})/offspring.total_units.^2;
    end
    
    % mutate activ fcns
    if config.multi_activ
        activFcn = offspring.activ_Fcn(i,:);
        pos =  randperm(length(activFcn),sum(rand(length(activFcn),1) < config.mut_rate));
        activFcn(pos) = {config.activ_list{randi([1 length(config.activ_list)],length(pos),1)}};
        offspring.activ_Fcn(i,:) = reshape(activFcn,size(offspring.activ_Fcn(i,:)));
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
        offspring.iir_weights{i,2}(pos,:) = [1+alpha -2*cos(w_0) 1-alpha];       
        
    end
end

% mutate output weights
if config.evolve_output_weights
    output_weights = offspring.output_weights(:);
    pos = randperm(length(output_weights),sum(rand(length(output_weights),1) < config.mut_rate));   
    output_weights(pos) = mutateWeight(output_weights(pos),[-config.output_weight_scaler config.output_weight_scaler],config);
    offspring.output_weights = reshape(output_weights,size(offspring.output_weights));
end

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