clear
close all

addpath(genpath('/mnt/lustre/users/md596/Magnets/MassMM')); % make sure filepath is correct

%set random seed for experiments
rng(1,'twister');

% type of network to evolve
config.res_type = 'MM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = [49];                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = 'narma_10';          % Task to evolve for
config.figure_array = [figure figure figure];
config.preprocess ='';

config.metrics = {'KR','GR','linearMC'};

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);
%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = 10;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = 1;

%% shape parameters
base_size = 125;
config.total_units = base_size*length(config.read_mag_direction);

shape_list = {'parallelogram','trapezoid','rectangle','custom','square','custom','circle'};
name_list = {'0','1','2','3','4','5','inf'};
shape_num_sides = [4, 4, 4, 3, 4, 5, 40];

theta = 0;
rot_point = [0.5 0.5];

% params sweep
data_points = 6;
IS_list = linspace(0.01,1,data_points); %logspace(1,3,5)/1e3
damping_list = linspace(0.1,1,data_points);%logspace(1,3,data_points)/1e3; %linspace(0.1,1,data_points); %

%cycle through shapes
for shape_indx = 1:length(shape_list)
    
    figure(config.figure_array(1))
    [x,y, num_points, area_ratio] = getPolyShape(shape_list{shape_indx}, shape_num_sides(shape_indx), theta, rot_point);
    total_cells_needed(shape_indx) = base_size/area_ratio;%(ceil(area_ratio*10)/10);%);
    square_film_dimensions = ceil(sqrt(total_cells_needed(shape_indx)));
    config.num_nodes = square_film_dimensions.^2;
    config.poly_num = num_points;
    
    % define individual
    default_pop = config.createFcn(config);
    
    % % reset input weights
    [xq,yq] = meshgrid(linspace(0,1,square_film_dimensions),linspace(0,1,square_film_dimensions));
    [in,on] = inpolygon(xq,yq,x,y);
    inputs_in_use = find(in | on);
    
    num_cells_used(shape_indx) = length(inputs_in_use);   
    
    for is_indx = 1:length(IS_list)
        
        for damping_indx = 1:length(damping_list)
            
            cnt = shape_indx*length(default_pop);
            
            parfor pop_indx = 1:length(default_pop)
                
                warning('off','all')
                rng(cnt + pop_indx,'twister');
                
                % assign coords
                default_pop(pop_indx).poly_coord = [x; y]';
                % apply inputs
                input_weights = sprandn(length(inputs_in_use),2, config.sparsity);
                default_pop(pop_indx).input_weights{1} = zeros(config.num_nodes,2);
                default_pop(pop_indx).input_weights{1}(inputs_in_use,:) = input_weights;               
                
                % params
                default_pop(pop_indx).input_scaling = IS_list(is_indx);
                default_pop(pop_indx).damping = damping_list(damping_indx);
                default_pop(pop_indx).leak_rate = 1;
                
                % evaluate
                %default_pop(pop_indx) = config.testFcn(default_pop(pop_indx),config);
                metrics(pop_indx,:) = getMetrics(default_pop(pop_indx),config);
                 
                fprintf('\n shape=%s , is=%.2f, damping=%.2f, i = %d, metrics = %d, %d, %.4f\n',shape_list{shape_indx},IS_list(is_indx),damping_list(damping_indx),pop_indx,metrics(pop_indx,:));
            end
              
            KR{shape_indx}(is_indx,damping_indx,:) =  metrics(:,1);
            GR{shape_indx}(is_indx,damping_indx,:) = metrics(:,2);
            MC{shape_indx}(is_indx,damping_indx,:) = metrics(:,3);

            avg_KR{shape_indx}(is_indx,damping_indx) = median(metrics(:,1));
            avg_GR{shape_indx}(is_indx,damping_indx) = median(metrics(:,2));
            avg_MC{shape_indx}(is_indx,damping_indx) = median(metrics(:,3));

            max_KR{shape_indx}(is_indx,damping_indx) = max(metrics(:,1));
            max_GR{shape_indx}(is_indx,damping_indx) = max(metrics(:,2));
            max_MC{shape_indx}(is_indx,damping_indx) = max(metrics(:,3));
   
            save('shape_metric_sweep.mat')
            
            
            % plot 
            h(1) = figure(config.figure_array(1));
            subplot(3,3,shape_indx)
            imagesc(avg_KR{shape_indx})
            colorbar
            title(shape_list(shape_indx))
            xticks(1:data_points)
            xticklabels(IS_list)
            xtickangle(45)
            yticks(1:data_points)
            yticklabels(damping_list)
            xlabel('IS')
            ylabel('Damping')
            caxis([0 200])
            drawnow
            
            h(2) = figure(config.figure_array(2));
            subplot(3,3,shape_indx)
            imagesc(avg_GR{shape_indx})
            colorbar
            title(shape_list(shape_indx))
            xticks(1:data_points)
            xticklabels(IS_list)
            xtickangle(45)
            yticks(1:data_points)
            yticklabels(damping_list)
            caxis([0 200])
            xlabel('IS')
            ylabel('Damping')
            drawnow
            
            h(3) = figure(config.figure_array(3));
            subplot(3,3,shape_indx)
            imagesc(avg_MC{shape_indx})
            title(shape_list(shape_indx))
            colorbar
            caxis([0 40])
             xticks(1:data_points)
            xticklabels(IS_list)
            xtickangle(45)
            yticks(1:data_points)
            yticklabels(damping_list)
            xlabel('IS')
            ylabel('Damping')
            drawnow
            
            savefig(h,'shape_metric.fig','compact')
        end
    end       
end