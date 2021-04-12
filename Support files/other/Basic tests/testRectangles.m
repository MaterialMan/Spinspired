clear
close all

%set random seed for experiments
rng(1,'twister');

%% Setup
config.parallel = 1;                        % use parallel toolbox

%start paralllel pool if empty
if isempty(gcp) && config.parallel %
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

% type of network to evolve
config.res_type = 'MM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = [100];                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = 'narma_10';          % Task to evolve for
config.figure_array = [figure figure];
config.preprocess ='';

config.metrics = {'KR','linearMC'}; 

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = 1;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = 1;

%% 
task_type = 0; % 1:tasks, 0: metrics
base_size = 100;
config.total_units = base_size*length(config.read_mag_direction);
shape_list = [2 4 5 10;...
              50 25 20 10]';
          
for shape = 1:size(shape_list,1)
    
    config.num_nodes = shape_list(shape,2)^2; % height
    default_pop = config.createFcn(config);

    mask = zeros(shape_list(shape,2),shape_list(shape,2));
    mask(1:shape_list(shape,1),1:shape_list(shape,2)) = 1; 
    inputs_in_use = find(mask);
    metrics = [];
    
    parfor input_pos = 1:length(inputs_in_use)
        warning('off','all')
        t_config = config;
        population(input_pos) = default_pop;
        population(input_pos).core_indx = input_pos;
        population(input_pos).input_scaling = 1;
        population(input_pos).damping = 0.1;
        population(input_pos).leak_rate = 1;
        population(input_pos).input_weights{1} = zeros(config.num_nodes,2);
        population(input_pos).input_weights{1}(inputs_in_use(input_pos),1) = 1;
        population(input_pos).geo_width = shape_list(shape,1)/shape_list(shape,2);
        population(input_pos).geo_height = 1;%shape_list(shape,1)/shape_list(shape,2);
        
        t = tic;
%         if task_type == 1
%             population(input_pos) = config.testFcn(population(input_pos),t_config);
%             error(input_pos) =  population(input_pos).test_error;
%             fprintf('Shape: %d, width: %.2f, height: %.2f, pos: %d, error: %.4f\n', shape,shape_list(shape,1),shape_list(shape,2),input_pos,error(input_pos))
%         else
            metrics(input_pos,:) = getMetrics(population(input_pos),t_config);
            fprintf('Shape: %d, width: %.2f, height: %.2f, pos: %d, Metrics: %.4f %.4f\n', shape,shape_list(shape,1),shape_list(shape,2),input_pos,metrics(input_pos,:))
        %end
       toc(t)
    end
    
    if task_type 
        figure(config.figure_array(1))
        set(gcf,'color','w');
        subplot(2,2,shape)
        error_matrix{shape} = zeros(size(mask));
        error_matrix{shape}(logical(mask(inputs_in_use))) = error;
        imagesc(error_matrix{shape})
        colorbar
        caxis([0 0.2])
        drawnow
        set(gca, 'ColorScale', 'log')
        save(strcat('shape_test_',config.dataset,'.mat'))
    else
        figure(config.figure_array(1))
        set(gcf,'color','w');
        subplot(2,2,shape)
        metrics_matrix{shape,1} = zeros(size(mask));
        metrics_matrix{shape,1}(logical(mask(inputs_in_use))) = metrics(:,1);
        imagesc(metrics_matrix{shape,1})
        colorbar
        caxis([0 config.total_units])
        
        figure(config.figure_array(2))
        set(gcf,'color','w');
        subplot(2,2,shape)
        metrics_matrix{shape,2} = zeros(size(mask));
        metrics_matrix{shape,2}(logical(mask(inputs_in_use))) = round(metrics(:,2));
        imagesc(metrics_matrix{shape,2})
        colorbar
        caxis([0 30])
        drawnow

        save(strcat('shape_test_metrics.mat'))
    end
    
end


figure
set(gcf,'color','w');
if task_type
    boxplot([error_matrix{1}(find(error_matrix{1})) error_matrix{2}(find(error_matrix{2})) error_matrix{3}(find(error_matrix{3})) error_matrix{4}(find(error_matrix{4}))],'Notch', 'on')
    set(gca, 'YScale', 'log')
    xlabel('Shape dimensions')
    ylabel('NMSE')
    xticklabels({'2 x 50', '4 x 25', '5 x 20', '10 x10' })
else
    subplot(1,2,1)
    boxplot([metrics_matrix{1,1}(find(metrics_matrix{1,1})) metrics_matrix{2,1}(find(metrics_matrix{2,1})) metrics_matrix{3,1}(find(metrics_matrix{3,1})) metrics_matrix{4,1}(find(metrics_matrix{4,1}))],'Notch', 'on')
    xlabel('Shape dimensions')
    ylabel('NMSE')
    xticklabels({'2 x 50', '4 x 25', '5 x 20', '10 x10' })

    subplot(1,2,2)
    boxplot([metrics_matrix{1,2}(logical(mask(inputs_in_use))) metrics_matrix{2,2}(find(metrics_matrix{2,2})) metrics_matrix{3,2}(find(metrics_matrix{3,2})) metrics_matrix{4,2}(find(metrics_matrix{4,2}))],'Notch', 'on')
    xlabel('Shape dimensions')
    ylabel('NMSE')
    xticklabels({'2 x 50', '4 x 25', '5 x 20', '10 x10' })
end


%% tot test
%states = config.assessFcn(population(input_pos),config.test_input_sequence,config,config.test_output_sequence);

%figure
% subplot(1,2,1)
% rectangle('Position',[0 0 1 shape_list(shape,1)/shape_list(shape,2)])
% axis([0 1 0 1])
% set(gca,'Ydir','reverse')
% 
% subplot(1,2,2)
% for t = 1:size(states,1)
%    x = reshape(states(t,config.num_nodes*2+1:end-1), shape_list(shape,2),shape_list(shape,2));
%    imagesc(x);
%    %set(gca,'Ydir','reverse')
%    caxis([-1 1]);
%    colormap(bluewhitered);
%    drawnow
% end
