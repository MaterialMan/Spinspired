function error = getError(error_to_check,individual)
switch(error_to_check)    
    case 'train'
        error = [individual.train_error];
    case 'val'
        error = [individual.val_error];
    case'test'
        error = [individual.test_error];
    case 'train&val'
        error = [individual.train_error] + [individual.val_error];
    case 'val&test'
        error = [individual.val_error] + [individual.test_error];
    case 'train&test'
        error = [individual.train_error] + [individual.test_error];
    case 'train&val&test'
        error = [individual.train_error] +[individual.val_error] +[individual.test_error];  
    case 'all&connectivity'
        for i = 1:length(individual)
             error(i) = individual(i).train_error + individual(i).val_error + individual(i).test_error...
                 + nnz(individual(i).input_weights)/(size(individual(i).input_weights,1)*size(individual(i).input_weights,2))...
                 + nnz(individual(i).restricted_output_mask)/size(individual(i).restricted_output_mask,1);%...
                % + nnz(individual(i).W)/(size(individual(i).W,1)*size(individual(i).W,2));
        end
    case 'test&connectivity'
         for i = 1:length(individual)
             error(i) = individual(i).test_error + nnz(individual(i).input_weights)/(size(individual(i).input_weights,1)*size(individual(i).input_weights,2));
         end
end