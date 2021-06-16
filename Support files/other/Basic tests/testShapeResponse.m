clear
close all

%set random seed for experiments
rng(1,'twister');

% type of network to evolve
config.res_type = 'MM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = [400];                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = 'test_sine';          % Task to evolve for
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

config.evolve_geometry = 1;                    % manipulate geomtry
config.evolve_poly = 1;
config.geometry_file =  'custom.geo';                    %add specific geometry file, e.g. custom.geo
config.lb = 0.1;
config.add_input_states = 0;

task_type = 0; % 1:tasks, 0: metrics
base_size = config.num_nodes;
config.total_units = base_size*length(config.read_mag_direction);
shape = 'circle';
num_sides = 4;
[x,y,num_points, area_ratio] = getPolyShape(shape,num_sides);

total_cells_needed = (base_size/area_ratio);
square_film_dimensions = round(sqrt(total_cells_needed));
config.num_nodes = square_film_dimensions.^2;
config.poly_num = num_points;

% define individual
default_pop = config.createFcn(config);
default_pop.poly_coord = [x; y]';
default_pop.core_indx = 1;
default_pop.input_scaling = 1;
default_pop.damping = 0.15;
default_pop.leak_rate = 1;

% % define inputs to sweep through
[xq,yq] = meshgrid(linspace(0,1,square_film_dimensions),linspace(0,1,square_film_dimensions));
[in,on] = inpolygon(xq,yq,x,y);
inputs_in_use = find(in | on);

default_pop.input_weights{1} = zeros(config.num_nodes,2);
switch(shape)
    %case {'rectangle'}
     %   input_pos = ceil(config.num_nodes/2) - ceil(sqrt(config.num_nodes)/4);
    %case {'rectangle','trapezoid','parallelogram'}
    %    input_pos = ceil(config.num_nodes/2) + ceil(sqrt(config.num_nodes)/4);
    %case 'triangle'
    %    input_pos = ceil(config.num_nodes/2)-1;
    otherwise
        if mod(sqrt(config.num_nodes),2) == 0
            input_pos = ceil(config.num_nodes/2) + ceil(sqrt(config.num_nodes)/2);
        else
            input_pos = ceil(config.num_nodes/2);
        end
end
default_pop.input_weights{1}(input_pos,1) = 1;

figure
s = zeros(square_film_dimensions);
s(inputs_in_use) = -1;
s(input_pos) = 1;
imagesc(s)
drawnow

states = config.assessFcn(default_pop,config.train_input_sequence,config);

figure
set(gcf,'color','w');
for t = 1:size(states,1)
    v = reshape(states(t,1:end), sqrt(config.num_nodes),sqrt(config.num_nodes));
    imagesc(v);
    %set(gca,'Ydir','reverse')
    caxis([min(min(states)) max(max(states))]);
    colormap(bluewhitered);
    %f(t) = getframe(gcf);
    drawnow
end

figure
set(gcf,'color','w');
subplot(2,2,1)
s = sum(states(:,1:end));
imagesc(reshape(s,sqrt(config.num_nodes),sqrt(config.num_nodes)))
colormap(gca,bluewhitered);
colorbar;
title('Sum')

subplot(2,2,2)
s = sum(abs(states(:,1:end)));
imagesc(reshape(s,sqrt(config.num_nodes),sqrt(config.num_nodes)))
colormap(gca,bluewhitered);
colorbar;
title('Sum (abs)')

subplot(2,2,3)
s = sum((states(:,1:end)).^2);
imagesc(reshape(s,sqrt(config.num_nodes),sqrt(config.num_nodes)))
colormap(gca,bluewhitered);
colorbar;
title('Sum squared')

subplot(2,2,4)
s = median(states(:,1:end));
imagesc(reshape(s,sqrt(config.num_nodes),sqrt(config.num_nodes)))
colormap(gca,bluewhitered);
colorbar;
title('Median')

figure
set(gcf,'color','w');
subplot(1,3,1)
s = sum(abs(states(:,1:config.num_nodes)));
imagesc(reshape(s,sqrt(config.num_nodes),sqrt(config.num_nodes)))
colormap(gca,bluewhitered);
colorbar;
title('X: Sum (abs)')

subplot(1,3,2)
s = sum(abs(states(:,1:config.num_nodes)));
imagesc(reshape(s,sqrt(config.num_nodes),sqrt(config.num_nodes)))
colormap(gca,bluewhitered);
colorbar;
title('Y: Sum (abs)')

subplot(1,3,3)
s = sum(abs(states(:,1:end)));
imagesc(reshape(s,sqrt(config.num_nodes),sqrt(config.num_nodes)))
colormap(gca,bluewhitered);
colorbar;
title('Z: Sum (abs)')

