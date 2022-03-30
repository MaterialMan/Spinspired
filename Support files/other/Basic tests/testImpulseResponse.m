clear
close all

%set random seed for experiments
rng(1,'twister');

%% Setup
config.parallel = 0;                        % use parallel toolbox

%start paralllel pool if empty
if isempty(gcp) && config.parallel %
    parpool('local',4,'IdleTimeout', Inf); % create parallel pool
end

% type of network to evolve
config.res_type = 'multiMM';                % state type of reservoir to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network of neurons), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more.
config.num_nodes = {[100]};                  % num of nodes in each sub-reservoir, e.g. if config.num_nodes = {10,5,15}, there would be 3 sub-reservoirs with 10, 5 and 15 nodes each. For one reservoir, sate as a non-cell, e.g. config.num_nodes = 25
config = selectReservoirType(config);   % collect function pointers for the selected reservoir type

%% Task parameters
config.discrete = 0;               % select '1' for binary input for discrete systems
config.nbits = 16;                 % only applied if config.discrete = 1; if wanting to convert data for binary/discrete systems
config.dataset = 'test_pulse';          % Task to evolve for
config.figure_array = [figure figure];
config.preprocess ='';

config.metrics = {'linearMC','KR','GR'}; 

% get any additional params stored in getDataSetInfo.m. This might include:
% details on reservoir structure, extra task variables, etc.
[config] = getAdditionalParameters(config);

% get dataset information
config = selectDataset(config);

%% Evolutionary parameters
config.num_tests = 1;                        % num of tests/runs
config.pop_size = 2;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)

%% sweep film sizes
config.test = 1;

% create population of reservoirs
population = config.createFcn(config);
default_pop = population(2);

%% damping effect
%damping = linspace(0.1,1,4);

time_scale = 0:10:200;
time_scale(1) = 1;

ppm = ParforProgMon('Initial population: ', length(time_scale));
for pop_indx = 1:1%length(time_scale)
    warning('off','all')

    tmp_config = config;
    population(pop_indx) = default_pop;
    %population(pop_indx).layer(1).input_scaling = 1;
    %population(pop_indx).layer(1).leak_rate = 1;
   % population(pop_indx).input_weights{1} = zeros(tmp_config.num_nodes,2);
   % population(pop_indx).input_weights{1}(tmp_config.num_nodes/2 - floor(sqrt(tmp_config.num_nodes)/2),1) = 1;
    
    population(pop_indx).layer(1).core_indx = pop_indx;
    population(pop_indx).layer(1).time_steps_increment = time_scale(pop_indx);
    
    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
    %MC(pop_indx,:) = getMetrics(population(pop_indx),tmp_config);
    
    %fprintf('Pop: %d, MC %.3f \n',pop_indx,MC(pop_indx))
    
    ppm.increment();
end

config.num_nodes = 100;

figure
%plot(time_scale,MC)
% set(gcf,'color','w')
for t = 1:size(states{1},1)
   imagesc(reshape(states{1}(t,1:end),sqrt(config.num_nodes),sqrt(config.num_nodes)));
   drawnow
end


d = 1;
[X,Y] = meshgrid(1:sqrt(config.num_nodes),1:sqrt(config.num_nodes));
for t = 1:size(states{d},1)
subplot(1,2,1)
surf(X,Y,reshape(states{d}(t,1:end-1),sqrt(config.num_nodes),sqrt(config.num_nodes))...
    ,'FaceAlpha',0.5,'FaceColor','interp');
set(gca,'visible','off')
%view([2 2 2])
zlim([min(min(states{d})) max(max(states{d}))])
caxis([-1 1])
xlabel('Side view')

subplot(1,2,2)
surf(X,Y,reshape(states{d}(t,1:end-1),sqrt(config.num_nodes),sqrt(config.num_nodes))...
    ,'FaceAlpha',0.5,'FaceColor','interp');
view(2)
caxis([-1 1])
set(gca,'visible','off')
xlabel('Top down view')
colorbar
colormap(bluewhitered)

drawnow;
end

figure
set(gcf,'color','w')
d = 1;
list = 11:2:15;
for t = 1:3
subplot(3,2,(t-1)*2+1)
surf(X,Y,reshape(states{d}(list(t),1:end-1),sqrt(config.num_nodes),sqrt(config.num_nodes))...
    ,'FaceAlpha',0.5,'FaceColor','interp');
set(gca,'visible','off')
%view([2 2 2])
zlim([-1 1])
caxis([-1 1])
if t ==1
    title('Side view')
end
xlabel(strcat('t = ',num2str(list(t))))
set(findall(gca,'type','text'),'visible','on')

subplot(3,2,t*2)
surf(X,Y,reshape(states{d}(list(t),1:end-1),sqrt(config.num_nodes),sqrt(config.num_nodes))...
    ,'FaceAlpha',0.5,'FaceColor','interp');
view(2)
caxis([-1 1])
set(gca,'visible','off')
if t ==1
    title('Top down view')
end
set(findall(gca,'type','text'),'visible','on')
colorbar
colormap(bluewhitered)
end

%% thickness vs time
colors = distinguishable_colors(length(time_scale));
figure
hold on
for i = 4:length(time_scale(4:end))
    p = plot(states{i}(:,201:300),'Color',colors(i,:));
    for n = 2:100
        set(get(get(p(n),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
hold off
xlabel('Time')
ylabel('Magnitude')
[~,hobj] = legend(string(round(time_scale(4:end)*100)/100));
h1 = findobj(hobj,'type','line');
set(h1,'LineWidth',4);

figure
bar(MC)
xticks(1:20)
xticklabels(round(time_scale*100)/100)
xtickangle(45)
xlabel('Damping')
ylabel('MC')

%% exchange effect
exchange = logspace(-26,-20,20);

ppm = ParforProgMon('Initial population: ', length(exchange));
parfor pop_indx = 1:length(exchange)
    warning('off','all')

    population(pop_indx) = default_pop;
    population(pop_indx).input_scaling = 1;
    population(pop_indx).leak_rate = 1;
    population(pop_indx).input_weights{1} = zeros(100,2);
    population(pop_indx).input_weights{1}(50,1) = 1;
    
    population(pop_indx).core_indx = pop_indx;
    population(pop_indx).damping = 0.5;
    population(pop_indx).exchange = exchange(pop_indx);
    
    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
    MC(pop_indx) = getMetrics(population(pop_indx),config);
    
    ppm.increment();
end

% thickness vs time
colors = distinguishable_colors(length(exchange));
figure
hold on
for i = 1:length(exchange)
    p = plot(states{i}(:,201:300),'Color',colors(i,:));
    for n = 2:length(p)
        set(get(get(p(n),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
hold off
xlabel('Time')
ylabel('Magnitude')
[~,hobj] = legend(string(exchange));
h1 = findobj(hobj,'type','line');
set(h1,'LineWidth',4);

figure
bar(MC)
xticks(1:20)
xticklabels(exchange)
xtickangle(45)
xlabel('Exchange')
ylabel('MC')

%% macro-cell size
cell_size = round(linspace(1,20,20)*10)/10;

ppm = ParforProgMon('Initial population: ', length(cell_size));
parfor pop_indx = 1:length(cell_size)
    warning('off','all')

    t_config = config;
    population(pop_indx) = default_pop;
    population(pop_indx).core_indx = pop_indx;
    
    population(pop_indx).input_scaling = 1;
    population(pop_indx).leak_rate = 1;
    population(pop_indx).system_size(1,1) = (sqrt(population(pop_indx).nodes(1)) * cell_size(pop_indx))-1;
    population(pop_indx).system_size(1,2) = (sqrt(population(pop_indx).nodes(1)) * cell_size(pop_indx))-1;
    t_config.macro_cell_size = cell_size(pop_indx);
    
    population(pop_indx).input_weights{1} = zeros(100,2);
    population(pop_indx).input_weights{1}(50,1) = 1;
    population(pop_indx).damping = 0.5;

    states{pop_indx} = config.assessFcn(population(pop_indx),t_config.test_input_sequence,t_config,t_config.test_output_sequence);
    MC(pop_indx) = getMetrics(population(pop_indx),t_config);
    
    ppm.increment();
end

% thickness vs time
colors = distinguishable_colors(length(cell_size));
figure
hold on
for i = 1:length(cell_size)
    p = plot(states{i}(:,201:300),'Color',colors(i,:));
    for n = 2:length(p)
        set(get(get(p(n),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
hold off
xlabel('Time')
ylabel('Magnitude')
[~,hobj] = legend(string(cell_size));
h1 = findobj(hobj,'type','line');
set(h1,'LineWidth',4);

figure
bar(MC)
xticks(1:20)
xticklabels(cell_size)
xtickangle(45)
xlabel('Macro-cell size')
ylabel('MC')

%% test exchange and damping
exchange = logspace(-26,-20,20);
time_scale = linspace(0,1,20);

ppm = ParforProgMon('Initial population: ', length(exchange)*length(time_scale));

for indx = 1:length(exchange)
    
    parfor pop_indx = 1:length(exchange)
        warning('off','all')
        
        population(pop_indx) = default_pop;
        population(pop_indx).input_scaling = 1;
        population(pop_indx).leak_rate = 1;
        population(pop_indx).input_weights{1} = zeros(100,2);
        population(pop_indx).input_weights{1}(50,1) = 1;
        
        population(pop_indx).core_indx = pop_indx;
        population(pop_indx).damping = time_scale(indx);
        population(pop_indx).exchange = exchange(pop_indx);
        
        %states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
        MC(indx,pop_indx) = getMetrics(population(pop_indx),config);
        
        ppm.increment();
    end
end



imagesc(MC)
xticks(1:20)
yticks(1:20)
xticklabels(exchange)
yticklabels(round(time_scale*100)/100)
xtickangle(45)
xlabel('Exchange')
ylabel('Damping')
a = colorbar;
a.Label.String = 'MC';

%% material density
density = round(linspace(0,1,20)*100)/100;

ppm = ParforProgMon('Initial population: ', length(density));
parfor pop_indx = 1:length(density)
    warning('off','all')

    population(pop_indx) = default_pop;
    population(pop_indx).input_scaling = 1;
    population(pop_indx).leak_rate = 1;
    population(pop_indx).input_weights{1} = zeros(100,2);
    population(pop_indx).input_weights{1}(50,1) = 1;
    
    population(pop_indx).core_indx = pop_indx;
    population(pop_indx).material_density = density(pop_indx);

    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
    MC(pop_indx) = getMetrics(population(pop_indx),config);
    
    ppm.increment();
end

% thickness vs time
colors = distinguishable_colors(length(density));
figure
hold on
for i = 2:length(density(1:end))
    p = plot(states{i}(:,201:300),'Color',colors(i,:));
    for n = 2:100
        set(get(get(p(n),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
hold off
xlabel('Time')
ylabel('Magnitude (x 1e-21)')
[~,hobj] = legend(string(density(2:end)));
h1 = findobj(hobj,'type','line');
set(h1,'LineWidth',4);

figure
bar(MC)
xticks(1:20)
xticklabels(density)
xtickangle(45)
xlabel('Density')
ylabel('MC')

%% applied field strength
field_strength = round(linspace(-1,1,20)*100)/100;
config.applied_field_unit_vector = {'0,0,1'};

ppm = ParforProgMon('Initial population: ', length(field_strength));
parfor pop_indx = 1:length(field_strength)
    warning('off','all')

    population(pop_indx) = default_pop;
    population(pop_indx).input_scaling = 1;
    population(pop_indx).leak_rate = 1;
    population(pop_indx).input_weights{1} = zeros(100,2);
    population(pop_indx).input_weights{1}(50,1) = 1;
    
    population(pop_indx).core_indx = pop_indx;
    population(pop_indx).applied_field_strength = field_strength(pop_indx);

    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
    MC(pop_indx) = getMetrics(population(pop_indx),config);
    
    ppm.increment();
end

% thickness vs time
colors = distinguishable_colors(length(field_strength));
figure
hold on
for i = 1:length(field_strength(1:end))
    p = plot(states{i}(:,201:300),'Color',colors(i,:));
    for n = 2:100
        set(get(get(p(n),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
hold off
xlabel('Time')
ylabel('Magnitude (x 1e-21)')
[~,hobj] = legend(string(field_strength(1:end)));
h1 = findobj(hobj,'type','line');
set(h1,'LineWidth',4);

figure
bar(MC)
xticks(1:20)
xticklabels(field_strength)
xtickangle(45)
xlabel('Applied field stength(T)')
ylabel('MC')

%% temperature
temperature = round(linspace(0,400,20));

ppm = ParforProgMon('Initial population: ', length(temperature));
parfor pop_indx = 1:length(temperature)
    warning('off','all')

    population(pop_indx) = default_pop;
    population(pop_indx).input_scaling = 1;
    population(pop_indx).leak_rate = 1;
    population(pop_indx).input_weights{1} = zeros(100,2);
    population(pop_indx).input_weights{1}(50,1) = 1;
    
    population(pop_indx).core_indx = pop_indx;
    population(pop_indx).temperature = temperature(pop_indx);

    states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
    MC(pop_indx) = getMetrics(population(pop_indx),config);
    
    ppm.increment();
end

% thickness vs time
colors = distinguishable_colors(length(temperature));
figure
hold on
for i = length(temperature(1:end)):-1:1
    p = plot(states{i}(:,201:300),'Color',colors(i,:));
    for n = 2:100
        set(get(get(p(n),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
hold off
xlabel('Time')
ylabel('Magnitude (x 1e-21)')
[~,hobj] = legend(string(temperature(end:-1:1)));
h1 = findobj(hobj,'type','line');
set(h1,'LineWidth',4);

figure
bar(MC)
xticks(1:20)
xticklabels(temperature)
xtickangle(45)
xlabel('Density')
ylabel('MC')

%% test temperature and film thickness
temperature = round(linspace(0,400,20));
thickness = round(linspace(0.1,2,10)*10)/10;

ppm = ParforProgMon('Initial population: ', length(temperature)*length(thickness));
MC = [];
for indx = 1:length(temperature)
    
    parfor pop_indx = 1:length(thickness)
        warning('off','all')
        
        population(pop_indx) = default_pop;
        population(pop_indx).input_scaling = 1;
        population(pop_indx).leak_rate = 1;
        population(pop_indx).input_weights{1} = zeros(100,2);
        population(pop_indx).input_weights{1}(50,1) = 1;
        
        population(pop_indx).core_indx = pop_indx;
        population(pop_indx).temperature = temperature(indx);
        population(pop_indx).system_size(3) = thickness(pop_indx);
        
        %states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
        MC(indx,pop_indx) = getMetrics(population(pop_indx),config);
        
        ppm.increment();
    end
end

figure
imagesc(MC')
xticks(1:length(temperature))
yticks(1:length(thickness))
xticklabels(temperature)
yticklabels(thickness)
xtickangle(45)
xlabel('Temperature (K)')
ylabel('Thickness (nm)')
a = colorbar;
a.Label.String = 'MC';

%% test temperature and macro size
temperature = round(linspace(0,400,20));
macro_size = round(linspace(1,15,10));

ppm = ParforProgMon('Initial population: ', length(temperature)*length(macro_size));
MC = [];
for indx = 1:length(temperature)
    
    parfor pop_indx = 1:length(macro_size)
        warning('off','all')
        
        t_config = config;
        population(pop_indx) = default_pop;
        population(pop_indx).input_scaling = 1;
        population(pop_indx).leak_rate = 1;
        population(pop_indx).input_weights{1} = zeros(100,2);
        population(pop_indx).input_weights{1}(50,1) = 1;
        
        population(pop_indx).core_indx = pop_indx;
        population(pop_indx).temperature = temperature(indx);
        
        population(pop_indx).system_size(1,1) = (sqrt(population(pop_indx).nodes(1)) * macro_size(pop_indx))-1;
        population(pop_indx).system_size(1,2) = (sqrt(population(pop_indx).nodes(1)) * macro_size(pop_indx))-1;
        t_config.macro_cell_size = macro_size(pop_indx);
        %states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
        MC(indx,pop_indx) = getMetrics(population(pop_indx),t_config);
        
        ppm.increment();
    end
end

figure
imagesc(MC')
xticks(1:length(temperature))
yticks(1:length(macro_size))
xticklabels(temperature)
yticklabels(macro_size)
xtickangle(45)
xlabel('Temperature (K)')
ylabel('Macro cell size (nm)')
a = colorbar;
a.Label.String = 'MC';

%% test temperature and input scaler
temperature = round(linspace(0,77,10));
input_scaler = [1 linspace(0,90,10)+10];

ppm = ParforProgMon('Initial population: ', length(temperature)*length(input_scaler));
MC = [];
for indx = 1:length(temperature)
    
    parfor pop_indx = 1:length(input_scaler)
        warning('off','all')
        
        t_config = config;
        population(pop_indx) = default_pop;
        population(pop_indx).input_scaling = 1*input_scaler(pop_indx);
        population(pop_indx).leak_rate = 1;
        population(pop_indx).input_weights{1} = zeros(100,2);
        population(pop_indx).input_weights{1}(50,1) = 1;
        
        population(pop_indx).core_indx = pop_indx;
        population(pop_indx).temperature = temperature(indx);
        

        states{pop_indx} = config.assessFcn(population(pop_indx),config.test_input_sequence,config,config.test_output_sequence);
        MC(indx,pop_indx) = getMetrics(population(pop_indx),t_config);
        
        ppm.increment();
    end
end

figure
imagesc(MC')
xticks(1:length(temperature))
yticks(1:length(input_scaler))
xticklabels(temperature)
yticklabels(input_scaler)
xtickangle(45)
xlabel('Temperature (K)')
ylabel('Input scaler')
a = colorbar;
a.Label.String = 'MC';