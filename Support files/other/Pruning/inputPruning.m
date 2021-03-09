% only works for one input currently
function [weights,rmv_cnt]= inputPruning(weights, config)

total_num_weights = size(weights,1)*size(weights,2);

total_num_connections_probability = nnz(weights)/total_num_weights;

input_matrix = reshape(full(weights(:,1)),sqrt(config.num_nodes),sqrt(config.num_nodes));
input_locations = find(input_matrix);
bias_matrix = reshape(full(weights(:,2)),sqrt(config.num_nodes),sqrt(config.num_nodes));
bias_locations = find(bias_matrix);

all_input_matrix = abs(input_matrix) + abs(bias_matrix);

% find all weights in use
in_use = find(all_input_matrix);
rmv_cnt = 0;    

for w = 1:length(in_use)
    
    % probability of pruning a weight based on number of weights in use
    if rand < total_num_connections_probability
        
        % calculate proxity to other weights
        [~,D] = knnsearch(in_use,in_use(w),'K',5);
        avg_dist = mean(D);
        proximity_probability = 1/avg_dist;
        
        if rand < proximity_probability
            % delete weight
            all_input_matrix(in_use(w)) = 0;  
            rmv_cnt = rmv_cnt+1;
        else            
            % calculate weight strength probability
            strength_probability = 1 - (abs(all_input_matrix(in_use(w)))/max(abs(all_input_matrix(in_use))));
            if rand < strength_probability
                % delete weight
                all_input_matrix(in_use(w)) = 0;
                rmv_cnt = rmv_cnt+1;
            end
           
        end
        
    end
    
end

% update weight matrices
input_matrix(input_locations) = all_input_matrix(input_locations);
bias_matrix(bias_locations) = all_input_matrix(bias_locations);

weights(:,1) = reshape(input_matrix,size(weights,1),1);
weights(:,2) = reshape(bias_matrix,size(weights,1),1);
    
    