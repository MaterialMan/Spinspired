function [current_individual,new_error]= pruneBlackbox(individual, config)

best_error = individual.test_error;
prev_individual = individual;

for iter = 1:config.prune_iterations
    
    current_individual = prev_individual;
    % prune weights
    [current_individual.input_weights,rmv_cnt]= inputPruning(current_individual.input_weights', config);
    current_individual.input_weights = current_individual.input_weights';
    %asses
    current_individual = config.testFcn(current_individual,config);
    % get new error
    new_error(iter) = current_individual.test_error;
        
    if best_error >= new_error(iter) && nnz(current_individual.input_weights) < nnz(prev_individual.input_weights)
        best_error = new_error(iter);
        prev_individual = current_individual;
        fprintf('Iter: %d, removed: %d, total Win: %d, last error: %.4f, new error: %.4f \n', iter, rmv_cnt, nnz(current_individual.input_weights), new_error(iter), best_error)
    end
end

end