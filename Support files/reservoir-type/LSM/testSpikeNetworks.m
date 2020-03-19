clear
load('tuning.mat')

% simulate the baseline period
[spikeMat_base, tVec_base] = poissonSpikeGen(6, 0.5, 20);
tVec_base = (tVec_base - tVec_base(end))*1000 - 1;
 
% simulate the stimulus period
[spikeMat_stim, tVec_stim] = poissonSpikeGen(30, 1, 20);
tVec_stim = tVec_stim*1000;
 
% put the baseline and stimulus periods together
spikeMat = [spikeMat_base spikeMat_stim];
tVec = [tVec_base tVec_stim];
 
% plot the raster and mark stimulus onset
plotRaster(spikeMat, tVec);
hold all;
plot([0 0], [0 size(spikeMat, 1)+1]);
 
% label the axes
xlabel('Time (ms)');
ylabel('Trial number');

% plot tuning curve of neuron 3
figure
hold on
for i = 1:size(tuningMat,1)
    plot(tuningMat(i,:))
end
hold off

%% Example 3 and 4
base_fr = 6;
stim_fr = 180;
num_trials = 30;

spikeMat_base =[];  spikeMat_stim =[];  spikeMat=[];

for i = 1:size(tuningMat,1)
    % simulate the baseline period
    [spikeMat_base{i}, tVec_base] = poissonSpikeGen(tuningMat(i,base_fr), 0.5, num_trials); 
    % simulate the stimulus period
    [spikeMat_stim{i}, tVec_stim] = poissonSpikeGen(tuningMat(i,stim_fr), 1, num_trials);  
    % put the baseline and stimulus periods together
    spikeMat{i} = logical([spikeMat_base{i} spikeMat_stim{i}]);
end

tVec_base = (tVec_base - tVec_base(end))*1000 - 1;
tVec_stim = tVec_stim*1000;
tVec = [tVec_base tVec_stim];
 
% plot the raster and mark stimulus onset
figure
for i = 1:size(tuningMat,1)
    subplot(4,2,i)
    plotRaster(spikeMat{i}, tVec);
    hold all;
    plot([0 0], [0 size(spikeMat{i}, 1)+1]);
end

 
