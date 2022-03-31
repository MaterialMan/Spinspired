function batch_path = getMMPath(batch_num)

%% check batch path - for multiple sessions of vampire using multiple cores (stops interference)
% add batch to path
current_batch = strcat('batch',num2str(batch_num),'.txt');
batch_path = which(current_batch);
batch_path = batch_path(1:end-length(current_batch)-1);

if isempty(batch_path) %&&  individual.core_indx == 1 % create new batch directory
    default_batch = 'batch0.txt'; % default batch
    default_batch_path = which(default_batch);
    path_to_copy = default_batch_path(1:end-length(default_batch)-1);
    new_dest_path = strcat(path_to_copy(1:end-length('batch0')),current_batch(1:end-4));
    
    % make new directory
    copyfile(path_to_copy,new_dest_path);
                    
    % add to path
    addpath(genpath(new_dest_path));
    
    % copy over files
    movefile(strcat(new_dest_path,'/',default_batch),strcat(new_dest_path,'/',current_batch))

    % assign new path
    batch_path = new_dest_path;
end

end