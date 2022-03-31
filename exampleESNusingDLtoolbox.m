
%% Japanese Vowels
clear all; close all; clc

[XTrain,TTrain] = japaneseVowelsTrainData;
inputSize = 12;
numHiddenUnits = 100;
numClasses = 9;

maxEpochs = 30;
miniBatchSize = 20;

% shorten padding of variable sequences
for i=1:numel(XTrain)
    sequence = XTrain{i};
    sequenceLengths(i) = size(sequence,2);
end

[sequenceLengths,idx] = sort(sequenceLengths,'descend');
XTrain = XTrain(idx);
TTrain = TTrain(idx);

% figure
% bar(sequenceLengths)
% xlabel("Sequence")
% ylabel("Length")
% title("Sorted Data")

acc = [];

for i = 1:1

    rng(i,"twister");
    % define network
    layers = [
        sequenceInputLayer(inputSize)
        RoRlayer(numHiddenUnits,inputSize,OutputMode="last")
        fullyConnectedLayer(numClasses)
        softmaxLayer
        classificationLayer];

    options = trainingOptions('adam',...
        'MaxEpochs',maxEpochs, ...
        'GradientThreshold',1, ...
        MiniBatchSize=miniBatchSize,...
        SequencePaddingDirection='left');
    %'Plots','training-progress',...

    %train
    net = trainNetwork(XTrain,TTrain,layers,options);

    % see outputs
    Y = predict(net,XTrain);

    % test
    [XTest,YTest] = japaneseVowelsTestData;
    YPred = classify(net,XTest,MiniBatchSize=miniBatchSize);
    %accuracy = mean(YTest==TTest)
    acc(i) = sum(YPred == YTest)./numel(YTest)

end

%% Chicken pox
clear all; close all; clc

data = chickenpox_dataset;
data = [data{:}];

% Partition data
numTimeStepsTrain = floor(0.8*numel(data));
dataTrain = data(1:numTimeStepsTrain+1);
dataTest = data(numTimeStepsTrain+1:end);

% rescale input
mu = mean(dataTrain);
sig = std(dataTrain);
dataTrainStandardized = (dataTrain - mu) / sig;
XTrain = dataTrainStandardized(1:end-1);
YTrain = dataTrainStandardized(2:end);

% Define network
numFeatures = 1;
numResponses = 1;
numHiddenUnits = 100;

layers = [ ...
    sequenceInputLayer(numFeatures)
    STOlayer(numHiddenUnits,numFeatures,OutputMode="sequence")
    fullyConnectedLayer(numResponses,'WeightsInitializer','he')
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',1, ...
    'GradientThreshold',2, ...
    'InitialLearnRate',0.1, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.2, ...
    'Plots','training-progress');

net = trainNetwork(XTrain,YTrain,layers,options);

dataTestStandardized = (dataTest - mu) / sig;
XTest = dataTestStandardized(1:end-1);

net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,YTrain(end));

numTimeStepsTest = numel(XTest);
for i = 2:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
end

YPred = sig*YPred + mu;

YTest = dataTest(2:end);
rmse = sqrt(mean((YPred-YTest).^2))

figure
plot(dataTrain(1:end-1))
hold on
idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
plot(idx,[data(numTimeStepsTrain) YPred],'.-')
hold off
xlabel("Month")
ylabel("Cases")
title("Forecast")
legend(["Observed" "Forecast"])

figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Forecast"])
ylabel("Cases")
title("Forecast")

subplot(2,1,2)
stem(YPred - YTest)
xlabel("Month")
ylabel("Error")
title("RMSE = " + rmse)