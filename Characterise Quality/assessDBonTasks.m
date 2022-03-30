%% Assess the database on tasks
% import behaviours and assess on tasks, output information as a input-output dataset 
function pred_dataset = assessDBonTasks(config,population,behaviours)

% Define tasks to evaluate
% Example: config.task_list = {'NARMA10','NARMA30','Laser','NonChanEqRodan'};
%parpool('local',4,'IdleTimeout', Inf); % create parallel pool

for set = 1:length(config.task_list)
        
    fprintf('Starting task: %s \n',config.task_list{set})
    % get datasets
    config.dataset = config.task_list{set};
    [config] = selectDataset(config);

    hw = waitbar(0,'Running...');
    %ppm = ParforProgMon('DB assessed: ', length(population));
    for indx = 1:length(population)
        population(indx) = config.testFcn(population(indx),config);
        train_error(indx,set) = population(indx).train_error;
        test_error(indx,set) = population(indx).test_error;

        if mod(indx,floor(length(population)/10))<1e-2
          waitbar(indx/length(population),hw);
        end
     
        %lineLength = fprintf('Indv %d, Train error: %.4f, test error: %.4f, on %s, ',indx,population(indx).train_error,population(indx).test_error,config.task_list{set});
       %pause(0.01)
        %fprintf(repmat('\b',1,lineLength))
      %  ppm.increment();
    end    
    hw =[];
end

% assign behaviours and task performance to struct for prediction
pred_dataset.inputs = behaviours;
pred_dataset.outputs.train_error = train_error;
pred_dataset.outputs.test_error = test_error;
       