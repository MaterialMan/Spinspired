% get STM of different sizes
function plotSTMofScalingSTOarrays(ID,num_points,end_grid_size)

addpath(genpath('/mnt/lustre/users/md596/Magnets/multilayer'));

% grid - change config.graph_type = {'basicLattice'}
array_type = {'basicLattice','Ring'};

% array types
for t = 1:length(array_type)
array_sizes = (6:2:end_grid_size).^2;
MC=[];

% all array sizes
for i = length(array_sizes):-1:1   
    close all
    population = STMFcn(ID,'STO',{array_sizes(i)},'MC',num_points,1,0,array_type{t});
    MC(i,:) = [population.metrics]; 
    close all
    save(strcat(array_type{t},'_STO_STMresults.mat'),'MC')
end
end

figure
boxplot(MC','notch','on')
