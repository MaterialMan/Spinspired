close all
rng(1,'twister')

figure1 =figure;
figure2=figure;
numSpinsPerDim = 100; %2^7
probSpinUp = 0.5;
%spin  = initSpins(numSpinsPerDim, probSpinUp);
J  = 0.5;
kT = 0;
cnt = 1;
%pos = numSpinsPerDim.^2/2 + numSpinsPerDim/2;
F = [];

% right = pos+repmat((numSpinsPerDim*(1:round(numSpinsPerDim*0.2))),round(numSpinsPerDim*0.2),1)+(1:round(numSpinsPerDim*0.2))';
% left = pos-repmat((numSpinsPerDim*(1:round(numSpinsPerDim*0.2))),round(numSpinsPerDim*0.2),1)-(1:round(numSpinsPerDim*0.2))';
% loc = [pos+(1:round(numSpinsPerDim*0.2)),pos-(1:round(numSpinsPerDim*0.2)), left(:)', right(:)'];
load('handDigits.mat');

for i = 1:100

%loc = randi([1 numSpinsPerDim.^2],round((numSpinsPerDim.^2)*0.75),1);
%spin  = initSpins(numSpinsPerDim, probSpinUp);
%spin(loc) = -spin(loc);

spin = ones(20);
spin((X(randi([1 length(X)]),:) > 0.5)) = -1;
spin = repelem(spin,5,5);

figure(figure1)
imagesc(spin);
drawnow

figure(figure2)
[spin,cnt]  = metropolis(spin, kT, J,cnt,probSpinUp);

fprintf('t = %d \n',i)
end
%Emean = energyIsing(spin, J);
%Mmean = magnetizationIsing(spin);