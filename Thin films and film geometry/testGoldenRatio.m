function testGoldenRatio(ID, task, pop_size)

close all

% print details
ID
task
pop_size

addpath(genpath('/mnt/lustre/users/md596/Magnets/MassMM')); % make sure filepath is correct

%set random seed for experiments
rng(ID,'twister');

%% Setup
config.parallel = 1;                        % use parallel toolbox

% type of network to evolve
config.res_type = 'MM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = [100];                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = task;          % Task to evolve for
config.preprocess ='';
config.figure_array = [figure];
config.metrics = {'KR','linearMC'};

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% Evolutionary parameters
config.num_tests = ID;                        % num of tests/runs
config.pop_size = pop_size;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = ID;

% define list of experiments
size_list = [8:2:16];

% test: 8x5 (40) ,10x6 (60), 12x7 (84), 14x9 (126), 16x10 (160)

%% INPUT MUST BE SYMMETRIC!!! only problem is this means sometimes more inputs and less outputs
Table = [];

%% Run experiments
for size_indx = 1:size(size_list,2)
    
    width  = size_list(size_indx);
    height = round(size_list(size_indx)/1.618);
    
    config.num_nodes = width.^2;
    
    default_pop = config.createFcn(config);
    
    % set dimensions
    tp_dim = height/width;
    dim_list = [1;tp_dim];
    
    metrics = [];
    error = [];
    
    % make rectangle
    coord = createRect(dim_list(2),dim_list(1));
    x = coord(:,1);
    y = coord(:,2);
    %area_ratio = polyarea(x',y');
    
    % % reset input weights
    square_film_dimensions = width;
    [xq,yq] = meshgrid(linspace(0,1,square_film_dimensions),linspace(0,1,square_film_dimensions));
    [in,on] = inpolygon(xq,yq,x,y);
    inputs_in_use = find(in | on);
    
    %define input
    if mod(square_film_dimensions,2) == 0
        %left (up & down)
        lu = round(length(inputs_in_use)/2)-floor(height/2);
        % make sure input is symmetrical in film
        if mod(height,2) == 0
            input_loc = inputs_in_use([lu lu+1 lu+height lu+height+1]); %[lu ld ru rd]
        else
            input_loc = inputs_in_use([lu lu+height]);
        end
    else
        lu = round(length(inputs_in_use)/2);
        % make sure input is symmetrical in film
        input_loc = inputs_in_use([lu lu+1]); %[lu ld ru rd]
    end
    
%     figure(config.figure_array(1))
%     subplot(2,2,4)
%     s = zeros(square_film_dimensions);
%     s(inputs_in_use) = -1;
%     s(input_loc) =1;
%     imagesc(s)
%     xticks(1:square_film_dimensions)
%     yticks(1:square_film_dimensions)
%     xticklabels(1:square_film_dimensions)
%     yticklabels(1:square_film_dimensions)
%     drawnow
    
    % stop broadcast variables
    geo_width = dim_list(2);
    geo_height = dim_list(1);
    size_in_use = width;
    size_div_in_use = height;
    
    parfor pop_indx = 1:config.pop_size
        warning('off','all')
        rng(pop_indx,'twister');
        
        t_config = config;
        %t_config.plot_states = 1;
        
        population(pop_indx) = default_pop(pop_indx);
        population(pop_indx).core_indx = pop_indx;
        
        population(pop_indx).input_scaling = 2*rand-1;
        population(pop_indx).damping = rand;
        population(pop_indx).leak_rate = rand;
        
        % input in centre of film
        population(pop_indx).input_weights{1} = zeros(population(pop_indx).nodes,  population(pop_indx).n_input_units + t_config.bias); %sprandn(population(pop_indx).nodes,  population(pop_indx).n_input_units + t_config.bias, t_config.sparsity);
        population(pop_indx).input_weights{1}(input_loc,1) = 1;
        
        population(pop_indx).geo_width = geo_width;
        population(pop_indx).geo_height = geo_height;
        
        t = tic;
        [population(pop_indx),test_states]= config.testFcn(population(pop_indx),t_config);
        % count no. of states used to check fair comparison
        num_states(pop_indx) = nnz(sum(test_states,1)) - config.add_input_states*population(pop_indx).n_input_units;
        population(pop_indx).total_units = num_states(pop_indx); %update to make KR measure manageable
        
        error(pop_indx) =  population(pop_indx).test_error;
        metrics(pop_indx,:) = getMetrics(population(pop_indx),t_config);
        population(pop_indx).behaviours = metrics(pop_indx,:);
        
        fprintf('Size: %d, width: %.2f, height: %.2f, pos: %d, error: %.4f, num cells: %d\n', size_in_use.^2,size_in_use,size_div_in_use,pop_indx,error(pop_indx),num_states(pop_indx))
        fprintf('.... Metrics: %.4f %.4f, took: %.3f\n',metrics(pop_indx,:),toc(t))
    end
    
    temp_table = [num_states;...                        % orignal size
        repmat(width,1,config.pop_size);...       % width
        repmat(height,1, config.pop_size);... % height
        [population.input_scaling];...
        [population.leak_rate];
        [population.damping];
        error;...                                                                           % error
        metrics']';                                                                          % metrics
    
    Table = [Table; temp_table];
    
    T = array2table(Table,'VariableNames',{'original_size','width','height','input_scaling','leak_rate','damping','error','KR','MC'});
    writetable(T,strcat('goldenRatio_',task,'_',num2str(pop_size),'_',num2str(ID),'.csv'),'Delimiter',',','QuoteStrings',true);
    
end

end