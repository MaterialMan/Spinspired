
%% Select Data Script: Generate task data sets and split data
function [config] = selectDataset(config)

scurr = rng;
temp_seed = scurr.Seed;

classification_data = 0;
wash_out = 50;

rng(1,'twister');
    
switch config.dataset
    
    %% test data
    case 'test_pulse'
        err_type = 'NMSE';
        wash_out =0;
        sequence_length = 300;
        config.train_fraction= 1/3;    config.val_fraction= 1/3;    config.test_fraction= 1/3;
        
        input_sequence = zeros(sequence_length,1);
        
        for i=1:sequence_length
            if mod(i,30) == 0
                input_sequence(i) = 1;
            end
        end
        
        ahead  = 1;
        output_sequence = input_sequence(1:end-ahead);
        input_sequence = input_sequence(ahead+1:end);
        
    case 'test_sequence'
        err_type = 'NMSE';
        wash_out =0;
        config.train_fraction= 1;    config.val_fraction= 0;    config.test_fraction= 0;
        
        input_sequence = [linspace(-10,10,1000) linspace(10,-10,1000)]';
        
        
        ahead  = 1;
        output_sequence = input_sequence(1:end-ahead);
        input_sequence = input_sequence(ahead+1:end);
        
        %% system modelling and chaotic systems
    case 'secondorder_task' %best 3.61e-3
        
        err_type = 'NMSE';
        
        sequence_length = 1500;
        config.train_fraction=0.5;    config.val_fraction=0.25;    config.test_fraction=0.25;
        
        u = rand(sequence_length,1)/2;
        y = zeros(sequence_length,1);
        for i = 3:sequence_length
            y(i) = 0.4*y(i-1)+0.4*y(i-1)*y(i-2)+0.6*(u(i).^2) + 0.1;
        end
        input_sequence = u;
        output_sequence = y;
        
    case 'narma_10' %input error 4 - good task
        err_type = 'NMSE';
        sequence_length = 5000;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,10);
        input_sequence = 2*input_sequence-0.5;
        output_sequence = 2*output_sequence-0.5;
        %fprintf('NARMA 10 task: %s \n',datestr(now, 'HH:MM:SS'))
        %config.preprocess = '';
        
    case 'narma_20' %input error 4 - good task
        err_type = 'NMSE';
        sequence_length = 5000;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,20);
        input_sequence = 2*input_sequence-0.5;
        output_sequence = 2*output_sequence-0.5;
        %fprintf('NARMA 20 task: %s \n',datestr(now, 'HH:MM:SS'))
        
        
    case 'narma_30' %input error 4 - good task
        err_type = 'NMSE';
        sequence_length = 5000;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,30);
        input_sequence = 2*input_sequence-0.5;
        output_sequence = 2*output_sequence-0.5;
        %fprintf('NARMA 30 task: %s \n',datestr(now, 'HH:MM:SS'))
        
    case 'multi_narma'
        err_type = 'NMSE';
        sequence_length = 5000;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        [input_sequence(:,1),output_sequence(:,1)] = generate_new_NARMA_sequence(sequence_length,5);
        [input_sequence(:,2),output_sequence(:,2)] = generate_new_NARMA_sequence(sequence_length,10);
        [input_sequence(:,3),output_sequence(:,3)] = generate_new_NARMA_sequence(sequence_length,30);
        
        input_sequence = 2*input_sequence-0.5;
        output_sequence = 2*output_sequence-0.5;
        
    case 'narma_10_DLexample' %input error 4 - good task
        err_type = 'NRMSE';
        config.preprocess = 0;
        sequence_length = 2700;
        config.train_fraction = 0.3333;    config.val_fraction=0.3333;    config.test_fraction=0.3333;
        
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,10);
        
    case 'narma_10_QRC' %input error 4 - good task
        err_type = 'NMSE';
        config.preprocess = 0;
        
        sequence_length = 12000;
        config.train_fraction = 0.3333;    config.val_fraction=0.3333;    config.test_fraction=0.3333;
        
        [input_sequence,output_sequence] = generate_new_NARMA_sequence(sequence_length,10);
        input_sequence = (input_sequence*2)*0.2;
        
    case 'henon_map' % input error > 1 - good task
%         
       % err_type = 'NMSE';
%         sequence_length= 8000;
%         stdev = 0.05;
%         config.train_fraction=0.375;    config.val_fraction=0.375;    config.test_fraction=0.25;
%         [input_sequence,output_sequence] = generateHenonMap(sequence_length,stdev);

        wash_out = 200;
        err_type = 'NRMSE_zhong';
        sequence_length= 2000;
        stdev = 0.05;
        config.train_fraction=0.5;    config.val_fraction=0;    config.test_fraction=0.5;
        [input_sequence,output_sequence] = generateHenonMap(sequence_length,stdev);
        
    case 'multi_signal'
        
        err_type = 'NMSE';
        
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        sequence_length = 2000;
        
        freq = 100;
        T = 100*(1/freq);
        fprintf('Freq: %d Hz\n',freq);
        config.Fs = 20000; %per channel
        step = 1/config.Fs;
        t = 0:step:T-step;
        amplitude = 1;
        t = t(:,1:sequence_length);
        
        % sinewave input
        input_sequence=[];
        input_sequence(:,1) = amplitude*sin(2*pi*freq*t);
        
        output_sequence = [];
        output_sequence(:,1) = amplitude*2*sawtooth(2*freq*t)-1;
        output_sequence(:,2) = amplitude*cos(2*pi*freq*t);
        output_sequence(:,3) = amplitude*square(2*pi*freq*t);
        output_sequence(:,4) = amplitude*sin(2*pi*2*freq*t);
        
        figure
        subplot(2,2,1)
        plot(input_sequence(1:1000,:))
        
        subplot(2,2,2)
        plot(output_sequence(1:1000,:))
        
        subplot(2,2,[3 4])
        xdft = fft(input_sequence);
        xdft = xdft(1:length(input_sequence)/2+1);
        freq = 0:config.Fs/length(input_sequence):config.Fs/2;
        plot(freq,abs(xdft));
        xlabel('Hz');
        
        %% Time-series
    case 'IPIX_plus5' % good task
        err_type = 'IPIX';
        sequence_length = 2000;
        config.train_fraction=0.4;    config.val_fraction=0.25;    config.test_fraction=0.35;   %val and test are switched later so ratios need to be swapped
        
        % IPIX radar task
        %load hiIPIX.txt
        load loIPIX.txt
        
        ahead = 5;
        data = loIPIX(1:sequence_length+ahead,:);
        input_sequence = data(1:sequence_length,:);
        output_sequence = data(ahead+1:end,:);
        
        %fprintf('Low IPIX task - 5 ahead. \n Started at %s \n',datestr(now, 'HH:MM:SS'))
        
    case 'IPIX_plus1' % good task
        err_type = 'IPIX';
        sequence_length = 2000;
        config.train_fraction=0.4;    config.val_fraction=0.25;    config.test_fraction=0.35;   %val and test are switched later so ratios need to be swapped
        
        % IPIX radar task
        %load hiIPIX.txt
        load loIPIX.txt
        
        ahead = 1;
        data = loIPIX(1:sequence_length+ahead,:);
        input_sequence = data(1:sequence_length,:);
        output_sequence = data(ahead+1:end,:);
        
        %fprintf('Low IPIX task 1 ahead. \n Started at %s \n',datestr(now, 'HH:MM:SS'))
        
        
    case 'laser' % good task
        
        err_type = 'NMSE';
        % Sante Fe Laser generator task
        sequence_length = 2000;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        ahead = 1;
        %data = laser_dataset;  %checkout the list at http://uk.mathworks.com/help/nnet/gs/neural-network-toolbox-sample-data-sets.html
        %data = cell2mat(data(:,1:sequence_length+ahead));
        data = load('laser.txt');
        input_sequence = data(1:sequence_length-ahead);
        output_sequence = data(ahead+1:sequence_length);
        
        %fprintf('Laser task TSP - 64 electrode test: %s \n',datestr(now, 'HH:MM:SS'))
        
    case 'sunspot' % good task but not sure about dataset- problem with dividing set
        
        % requires deep learning toolbox
        err_type = 'NMSE';
        % Sunspot task - needs proper dataset separation
        sequence_length = 2899;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        ahead = 1;
        %load sunspot.txt %solar_dataset;  %checkout the list at http://uk.mathworks.com/help/nnet/gs/neural-network-toolbox-sample-data-sets.html
        data = solar_dataset;
        data = [data{:}]';
        
        %data = sunspot(1:sequence_length+ahead,4);
        input_sequence = data(1:sequence_length-ahead);
        output_sequence = data(ahead+1:sequence_length);
        
        %fprintf('Sunspot task TSP: %s \n',datestr(now, 'HH:MM:SS'))
        
        %% Learning attractors
    case 'spike_attractor'
        
        config.preprocess = '';
        config.preprocess_shift = [0 1];
        
        err_type = 'MSE';
        
        wash_out =50;
        sequence_length = 4000+1000;
        config.train_fraction=0.8;    config.val_fraction=0;    config.test_fraction=0.2;
        
        %multiple attractor dataset
        num_attractors = 20;
        
        % train sequence
        train_u = zeros(sequence_length/2,num_attractors);
        train_y = ones(sequence_length/2,num_attractors)*-0.5;
        
        train_u(1,1) = 0.5;
        indx=1;
        for i = 1:length(train_u)
            if mod(i,200) ==0
                
                if indx == num_attractors
                    indx=1;
                else
                    indx = indx +1;
                end
                train_u(i+1,indx) = 0.5;
            end
            train_y(i,indx) = 0.5;
        end
        
        % test sequence
        test_u = zeros(sequence_length/2,num_attractors);
        test_y = ones(sequence_length/2,num_attractors)*-0.5;
        
        indx = [];
        for i = 1:length(test_u)
            if rand < 0.02
                indx = randi([1 num_attractors]);
                test_u(i,indx) = 0.5;
            end
            
            if ~isempty(indx)
                test_y(i,indx) = 0.5;
            end
        end
        
        input_sequence = [train_u; test_u];
        output_sequence = [train_y; test_y];
        
    case 'lorenz'
        err_type = 'NMSE';
        config.train_fraction=1;    config.val_fraction=0;    config.test_fraction=0;
        wash_out =0;
        
        data_length = 1e3; T = 100; h = 0.001;
        [x,y, z] = createLorenz(28, 10, 8/3, T, h, data_length);  % roughly 100k datapoints
        input_sequence= [1,1,1 ;zeros(data_length-1,3)];
        output_sequence= [x, y, z];
        
    case 'attractor' %reconstruct lorenz attractor
        err_type = 'NMSE';
        config.train_fraction=0.3333;    config.val_fraction=0;    config.test_fraction=0.6666;
        wash_out =100;
        
        switch(config.attractor_type)
            case 'lorenz'
                data_length = 4e3; T = 100; h = 0.001;
                [x,y, z] = createLorenz(28, 10, 8/3, T, h, data_length);  % roughly 100k datapoints
                input_sequence= [x, y, z];
                slice = [1 1 0.5];
            case 'rossler'
                data_length = 4e3; T = 100; h = 0.001;
                [x,y,z] = createRosslerAttractor(0.2,0.2,5.7, T, h ,data_length); % roughly 100k datapoints
                input_sequence= [x, y, z];
                slice = [1 1 0.5];
            case 'limit_cycle'
                data_length = 4e3; T = 100; h = 0.001;
                [x, y] = createLimitCycleAttractor(4, T, h, data_length); % roughly 10k datapoints
                input_sequence= [x, y];
                slice = [1 1 0.5];
            case 'mackey_glass'
                data_length = 3e3; T = 1e4;
                [x] = createMackeyGlass(17, 0.1, 0.2, 10, T ,data_length);
                input_sequence= x';
                %x = load('Mackey_Glass_t17.txt');
                %attractor_sequence= x(1:data_length);
                slice = [1 1 0.5];
            case 'duffing_map'
                data_length = 4e3;
                data_struct.delta= 0.3;
                data_struct.alpha= -1;
                data_struct.beta= 1;
                data_struct.gamma= 0.5;
                data_struct.w = 1.2;
                y0 = [1 0];
                T = 1e3;
                [x] =createDuffingOscillator(data_length, data_struct, y0, T);
                input_sequence= x';
                slice = [1 1 0.5];
            case 'dynamic' % not finished: still playing with
                data_length = 4e3;
                plot_on = 1;
                num_attractors = 10;
                [x] = attractorSwitch(data_length,num_attractors,plot_on);
                input_sequence= x;
                slice = [1 1 1];
            case 'multi_attractor'
                err_type = 'NMSE';
                %multiple attractor dataset
                sequence_length = 5000; %4000 train, 1000 test
                num_attractors = 3;
                config.train_fraction=0.8;    config.val_fraction=0;    config.test_fraction=0.2;
                
                % create train sequence
                train_u = zeros(sequence_length*config.train_fraction,num_attractors);
                train_y = ones(sequence_length*config.train_fraction,num_attractors)*-0.5;
                
                train_u(1,1) = 0.5;
                indx=1;
                for i = 1:length(train_u)
                    if mod(i,200) ==0 && i<length(train_u)
                        if indx == num_attractors
                            indx=1;
                        else
                            indx = indx +1;
                        end
                        train_u(i+1,indx) = 0.5;
                    end
                    train_y(i,indx) = 0.5;
                end
                
                % test sequence
                test_u = zeros(sequence_length*config.test_fraction,num_attractors);
                test_y = ones(sequence_length*config.test_fraction,num_attractors)*-0.5;
                
                indx = [];
                for i = 1:length(test_u)
                    if rand < 0.02
                        indx = randi([1 num_attractors]);
                        test_u(i,indx) = 0.5;
                    end
                    
                    if ~isempty(indx)
                        test_y(i,indx) = 0.5;
                    end
                end
                
                input_sequence = [train_u; test_u];
                output_sequence = [train_y; test_y];
                
                slice = [1 1 1];
            otherwise
        end
        
        data_length = size(input_sequence,1);
        
        if ~strcmp(config.attractor_type,'multi_attractor')
            ahead = 1;%shift by 1. Becomes prediction problem
            base_sequence = input_sequence;
            input_sequence = base_sequence(1:end-ahead,:);
            output_sequence = base_sequence(1+ahead:end,:);
        end
        
        % divide data -  add no signal
        input_sequence(floor(data_length*config.train_fraction*slice(1))+1:floor(data_length*config.train_fraction),:) = zeros;
        input_sequence(floor(data_length*config.train_fraction)+floor(data_length*config.val_fraction*slice(2))+1:floor(data_length*config.train_fraction)+floor(data_length*config.val_fraction),:) = zeros;
        input_sequence(floor(data_length*config.train_fraction)+floor(data_length*config.val_fraction)+floor(data_length*config.test_fraction*slice(3))+1:end,:) = zeros;
        
        %% Pattern Recognition - using PCA to reduce dimensions maybe very useful
    case 'MNIST'
        
        err_type = 'softmax';
        config.preprocess = 0;
        config.preprocess_shift = 'zero to one';
        data_length = 2000;
        
        wash_out = 0;
        config.train_fraction=0.8;    config.val_fraction=0.1;    config.test_fraction=0.1;
        
        [XTrain, YTrain, XTest, YTest] = load_mnist('./Support files/other/Datasets/MNIST');
        
        input_sequence = [XTrain; XTest];
        output_sequence = [YTrain; YTest];
        
        % make sure distribution is even
        indx = [];
        for k = 1:10
            n = find(output_sequence(:,k));
            indx(:,k) = n(randi([1 length(n)],1,data_length/10));
        end
        % get correct lengths and id's for data split
        train_id = indx(1:(data_length/10)*config.train_fraction,:);
        val_id = indx(1:(data_length/10)*config.val_fraction,:);
        test_id = indx(1:(data_length/10)*config.test_fraction,:);
        
        % randomise
        train_id = train_id(randperm(length(train_id(:))));
        val_id = val_id(randperm(length(val_id(:))));
        test_id = test_id(randperm(length(test_id(:))));
        
        input_sequence = input_sequence([train_id val_id test_id],:);
        output_sequence = output_sequence([train_id val_id test_id],:);
        
    case 'NIST-64' %Paper: Reservoir-based techniques for speech recognition
        err_type = 'OneVsAll_NIST';
        xvalDetails.kfold = 5;
        xvalDetails.kfoldSize = 150;
        xvalDetails.kfoldType = 'standard';
        config.train_fraction=0.7;    config.val_fraction=0.15;    config.test_fraction=0.15; %meaningless
        
        y_list = [];
        u_list = [];
        lens = [];
        
        for i = 1:5
            l = [ 1 2 5 6 7];
            for j = 1:10
                for n = 0:9
                    u_z = zeros(77,xvalDetails.kfoldSize);
                    u = load(strcat('s',num2str(l(i)),'_u',num2str(j),'_d',num2str(n)));
                    u_z(:,1:size(u.spec,2)) = u.spec;
                    
                    y = zeros(10,size(u_z,2))-1;
                    y(n+1,:) = ones(1,size(u_z,2));
                    u_list = [u_list u_z];
                    lens = [lens size(u_z,2)];
                    y_list = [y_list y];
                end
            end
        end
        
        input_sequence = u_list';
        output_sequence = y_list';
        
    case 'hand_digits'
        
        wash_out = 0;
        err_type = 'softmax';
        config.preprocess = 0;
        config.preprocess_shift = 'zero to one';
        config.train_fraction=0.8;    config.val_fraction=0.1;    config.test_fraction=0.1;
        
        dataset_length = 5000; %manually change dataset length for xval
        
        load('handDigits.mat');
        input_sequence = X;
        output_sequence = [];
        for i = 1:10
            output_sequence(:,i) = y==i;
        end
        
        target=randperm(dataset_length);
        temp_inputSequence = input_sequence(target,:);
        temp_outputSequence = output_sequence(target,:);
        
        input_sequence = temp_inputSequence;
        output_sequence = temp_outputSequence;
        
    case 'cifar10'
        
        err_type = 'softmax';
        config.preprocess = 0;
        config.preprocess_shift = 'zero to one';
        data_length = 5000;
        
        wash_out = 0;
        config.train_fraction=0.8;    config.val_fraction=0.1;    config.test_fraction=0.1;
        
        [XTrain, YTrain, XTest, YTest] = load_cifar10;
        
        input_sequence = [XTrain; XTest];
        output_sequence = [YTrain; YTest];
        
        % make sure distribution is even
        indx = [];
        for k = 1:10
            n = find(output_sequence(:,k));
            indx(:,k) = n(randi([1 length(n)],1,data_length/10));
        end
        % get correct lengths and id's for data split
        train_id = indx(1:(data_length/10)*config.train_fraction,:);
        val_id = indx(1:(data_length/10)*config.val_fraction,:);
        test_id = indx(1:(data_length/10)*config.test_fraction,:);
        
        % randomise
        train_id = train_id(randperm(length(train_id(:))));
        val_id = val_id(randperm(length(val_id(:))));
        test_id = test_id(randperm(length(test_id(:))));
        
        input_sequence = input_sequence([train_id val_id test_id],:);
        output_sequence = output_sequence([train_id val_id test_id],:);
        
    case 'japanese_vowels' %(12: IN, 9:OUT - binary ) - input only 83% accuracy!  Train:0.2288  Test:0.1863
        err_type = 'softmax'; %Paper: Optimization and applications of echo state networks with leaky- integrator neurons
        
        % Nine male speakers uttered two Japanese vowels /ae/ successively.
        % For each utterance, with the analysis parameters described below, we applied
        % 12-degree linear prediction analysis to it to obtain a discrete-time series
        % with 12 LPC cepstrum coefficients. This means that one utterance by a speaker
        % forms a time series whose length is in the range 7-29 and each point of a time
        % series is of 12 features (12 coefficients).
        % The number of the time series is 640 in total. We used one set of 270 time series for
        % training and the other set of 370 time series for testing.
        
        [train_input_sequence,trainOutputSequence,testInputSequence,test_output_sequence] = readJapVowels();
        input_sequence = [train_input_sequence; testInputSequence];
        output_sequence = [trainOutputSequence; test_output_sequence];
        config.train_fraction=size(train_input_sequence,1)/9961;    config.val_fraction=(size(testInputSequence,1)/9961)*0.1;    config.test_fraction=(size(testInputSequence,1)/9961)*0.9;
        
        %        t =  randperm(dataset_length,dataset_length);
        
        %% signal recovery
    case 'non_chan_eq_rodan' % (1:in, 1:out) error 0.999 Good task, requires memory
        err_type = 'NMSE';
        %input alone error = 0.091
        sequence_length = 2000;
        config.train_fraction=0.25;    config.val_fraction=0.375;    config.test_fraction=0.375;
        
        [input_sequence, output_sequence] = NonLinear_ChanEQ_data(sequence_length);
        input_sequence =input_sequence';
        output_sequence =output_sequence';
        
        %% classification
    case 'signal_classification'
        err_type = 'softmax';
        
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        freq = 1000;
        fprintf('Signal Classification: \n',datestr(now, 'HH:MM:SS'))
        fprintf('Freq: %d Hz\n',freq);
        scanFreq = 20000; %per channel
        step = 1/scanFreq;
        t = 0:step:1-step;
        amplitude = 1;
        sequence_length = 3000;
        period = 20;
        
        % sinewave input
        input_sequence(:,1) = amplitude*sin(2*pi*freq*t);
        input_sequence(:,2) = amplitude*square(2*pi*freq*t);
        
        cnt = 1; sinInput =[];squareInput=[];
        for i = 0:period:sequence_length-period
            sinInput(cnt,i+1:i+period) = input_sequence(i+1:i+period,1);
            squareInput(cnt,i+1:i+period) = input_sequence(i+1:i+period,2);
            cnt = cnt +1;
        end
        
        combInput = zeros(sequence_length,1);
        combOutput= ones(sequence_length,2)*0;
        for i = 1:sequence_length/period
            if round(rand)
                combInput = combInput+sinInput(i,:)';
                combOutput((i*period)-period+1:i*period,1) =  ones(period,1);
            else
                combInput = combInput+squareInput(i,:)';
                combOutput((i*period)-period+1:i*period,2) =  ones(period,1);
            end
        end
        
        input_sequence = combInput;
        output_sequence = combOutput;
        
    case 'spiral'
        classification_data = 1;
        %config.preprocess = 0;
        %config.preprocess_shift = 'zero to one';
        
        err_type = 'softmax';
        
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        num_samples = 400;
        noise = 0.1;
        
        [input_sequence, output] = createSpiral(num_samples,noise);
        
        output_sequence = zeros(length(input_sequence),2);
        output_sequence(output==-1,1) = 1;
        output_sequence(output==1,2) = 1;
        
        t =  randperm(length(input_sequence),length(input_sequence));
        
        input_sequence = input_sequence(t,:);
        output_sequence = output_sequence(t,:);
        
        scatter(input_sequence(:,1), input_sequence(:,2), 15, output_sequence(:,1), 'filled')
        
    case 'spiral_multi_class'
        
        classification_data = 1;
        %config.preprocess = 0;
        %config.preprocess_shift = 'zero to one';
        
        err_type = 'softmax';
        wash_out = 0;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        N = 400; % number of points per class
        D = 2; % dimensionality
        K = 5; % number of classes
        X = zeros(N*K,D); % data matrix (each row = single example)
        y = zeros(N*K,K); % class labels
        for j = 0:K-1
            ix = N*j+1:N*(j+1);
            r = linspace(0,1,N); % radius
            t = linspace(j*(5),(j+1)*(5),N) + randn(1,N)*0.2; % theta
            X(ix,:) = [r.*sin(t); r.*cos(t)]';
            y(ix,j+1) = 1;
        end
        % lets visualize the data:
        scatter(X(:,1), X(:,2), 40)
        
        t =  randperm(length(X),length(X));
        
        input_sequence = X(t,:);
        output_sequence = y(t,:);
        
        
    case 'iris' %iris_dataset; (4:in, 3:out) %input alone 76% - medium task
        
        classification_data = 1;
        
        err_type = 'softmax';%'IJCNNpaper';%'confusion';
        
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        dataset_length = 150;
        
        t =  randperm(dataset_length,dataset_length);
        
        load('iris.mat');
        
        %[input_sequence, output_sequence] =  %iris_dataset; %iris_dataset; (4:in, 3:out)
        input_sequence = input_sequence(:,t)';
        output_sequence = output_sequence(:,t)';
        
    case 'breast_cancer'
        
        classification_data = 1;
        % breast cancer
        err_type = 'softmax';
        
        
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        %[input_sequence, output_sequence] =  cancer_dataset;
        
        T = readtable('wisconsin_WDBC.csv', 'HeaderLines',1);
        T = T{:,:};
        
        input_sequence = T(:,1:end-1);
        output=T(:,end);
        output_sequence = zeros(length(input_sequence),2);
        output_sequence(output==2,2) = 1;
        output_sequence(output==1,1) = 1;
        
        t =  randperm(size(input_sequence,1),size(input_sequence,1));
        input_sequence = input_sequence(t,:);
        output_sequence = output_sequence(t,:);
        
    case 'robot-wall-follow'
        
        classification_data = 1;
        % breast cancer
        err_type = 'softmax';
        
        
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        %[input_sequence, output_sequence] =  cancer_dataset;
        
        T = readtable('robot-wall-follow-classification.csv', 'HeaderLines',1);
        T = T{:,:};
        
        input_sequence = T(:,1:end-1);
        output=T(:,end);
        output_sequence = zeros(length(input_sequence),4);
        output_sequence(output==1,1) = 1;
        output_sequence(output==2,2) = 1;
        output_sequence(output==3,3) = 1;
        output_sequence(output==4,4) = 1;
        
        t =  randperm(size(input_sequence,1),size(input_sequence,1));
        input_sequence = input_sequence(t,:);
        output_sequence = output_sequence(t,:);
        
        
        
        %% Multi-timescale
    case {'MSO1','MSO2','MSO3','MSO4','MSO5','MSO6','MSO7','MSO8','MSO9','MSO10','MSO11','MSO12'}  %MSO'
        
        %config.preprocess = '';
        %config.preprocess_shift = [0 1]; % range for data
        
        task = str2num(config.dataset(4:end));
        err_type = 'NRMSE';
        wash_out = 100;
        sequence_length= 1000;
        config.train_fraction=0.4;    config.val_fraction=0.3;    config.test_fraction=0.3;
        
        ahead = 1;
        for t = 1:sequence_length+ahead
            u(t,1) = sin(0.2*t);
            u(t,2) = u(t,1) + sin(0.311*t);
            u(t,3) = u(t,2) + sin(0.42*t);
            u(t,4) = u(t,3) + sin(0.51*t);
            u(t,5) = u(t,4) + sin(0.63*t);
            u(t,6) = u(t,5) + sin(0.74*t);
            u(t,7) = u(t,6) + sin(0.85*t);
            u(t,8) = u(t,7) + sin(0.97*t);
            u(t,9) = u(t,8) + sin(1.08*t);
            u(t,10) = u(t,9) + sin(1.19*t);
            u(t,11) = u(t,10) + sin(1.27*t);
            u(t,12) = u(t,11) + sin(1.32*t);
        end
        
        % freq = [5 50 500 5000 50000];
        Fs = 2000;            % Sampling frequency
        % T = 1/Fs;             % Sampling period
        % L = 2000;             % Length of signal
        % t = (0:L-1)*T;        % Time vector
        % u(:,1) =  sin(2*pi*freq(1)*t)';
        % for n = 2:5
        %     u(:,n) = u(:,n-1) + sin(2*pi*freq(n)*t)';
        % end
        
        L = sequence_length;
        subplot(1,2,1)
        xdft = fft(u(:,task));
        P2 = abs(xdft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        plot(f,P1)
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        subplot(1,2,2)
        plot(u(:,task))
        
        % remember
        % input_sequence = u(ahead+1:sequence_length+ahead,task);
        % output_sequence = u(1:sequence_length,task);
        
        % predict
        input_sequence = u(1:sequence_length,task);
        output_sequence = u(ahead+1:sequence_length+ahead,task);
        
        % [output_sequence] = featureNormailse(output_sequence,config);
        
        %% Reinforcement - no teacher signal
    case 'pole_balance'
        err_type = 'empty';
        config.train_fraction=0.1;    config.val_fraction=0.1;    config.test_fraction=0.1;
        wash_out = 1;
        input_sequence= zeros(100,6);
        output_sequence= zeros(100,1);
        
        
    case 'robot'
        err_type = 'empty';
        config.train_fraction=0.1;    config.val_fraction=0.1;    config.test_fraction=0.1;
        wash_out = 1;
        config.num_sensors = 8;
        input_sequence= zeros(100,config.num_sensors+1);
        output_sequence= zeros(100,4);
        
        
        %% Binary, adders, logic etc
    case 'binary_nbit_adder'
        
        err_type = 'hamming';
        type = 'nbit_adder';
        bit = 3;
        datalength = 5000;
        config.train_fraction=0.5;    config.val_fraction=0.25;    config.test_fraction=0.25;
        config.preprocess =0;
        
        A_in = randi([0 (2^bit)-1],datalength,1);
        B_in = randi([0 (2^bit)-1],datalength,1);
        
        input = [de2bi(A_in,bit) de2bi(B_in,bit)];
        output=[];
        for i = 1:datalength
            output(i,:) = getAdderTruthTable(type,[bit,A_in(i),B_in(i)]);
        end
        
        in  = [bi2de(input(:,1:bit)) bi2de(input(:,bit+1:end))];
        out = bi2de(output);
        hist(out)
        
        % get uniform distribution
        [N,edges,bin] = histcounts(out,2^bit*2);
        cnt = 1;bin_in=[];bin_out=[];
        for i = min(bin):max(bin)
            bin_in{cnt} = input(bin == i,:);
            bin_out{cnt} = output(bin == i,:);
            cnt = cnt +1;
        end
        
        for i = 1:datalength
            ex = 1;
            while(ex)
                pos = randi([min(bin) max(bin)]);
                if ~isempty(bin_in{pos})
                    ex = 0;
                end
            end
            pos2 = randi([1 length(bin_in{pos})]);
            input_sequence(i,:) = bin_in{pos}(pos2,:);
            output_sequence(i,:) = bin_out{pos}(pos2,:);
        end
        
        in  = [bi2de(input_sequence(:,1:bit)) bi2de(input_sequence(:,bit+1:end))];
        out = bi2de(output_sequence);
        hist(out)
        %hist(in(:,2))
        %hist(in(:,1))
        
        %% Filtering
    case 'image_gaussian' % Gaussian noise task
        err_type = 'NMSE';
        wash_out = 0;
        config.train_fraction= 0.5;    config.val_fraction=0.25;    config.test_fraction=0.25;
        image_size = 32;
        grey = 1;
        
        % load('airplanes_800x25x25.mat');
        input_sequence = []; output_sequence =[];
        for i = 1:800
            if i > 9 && i < 100
                file_dir = strcat('airplanes\image_00',num2str(i),'.jpg');
            elseif i > 99
                file_dir = strcat('airplanes\image_0',num2str(i),'.jpg');
            else
                file_dir = strcat('airplanes\image_000',num2str(i),'.jpg');
            end
            [img1_input,img1_output] = getImage(file_dir, image_size, 'gaussian',grey);
            
            input_sequence = [input_sequence; img1_input(:)'];
            output_sequence = [output_sequence; img1_output(:)'];
        end
        
        %% Interpolation and function fitting
    case 'franke_fcn' % Task is to emulate a nonlinear function
        err_type = 'NMSE';
        sequence_length = 1024;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        data = linspace(0,1,sqrt(sequence_length));
        [X,Y] = meshgrid(data);
        input_sequence = [X(:),Y(:)];
        
        
        for i = 1:sequence_length
            output_sequence(i) = franke2d(input_sequence(i,1),input_sequence(i,2));
        end
        
        t =  randperm(sequence_length);
        input_sequence = input_sequence(t,:);
        output_sequence = output_sequence(t)';
        
        scatter3(input_sequence(:,1),input_sequence(:,2),output_sequence)
        
    case 'chemical_fit'
        %     8x498 matrix: measurements taken from eight sensors during a chemical process.
        %     1x498 matrix of a ninth sensor's measurements, to be estimated from the first eight.
        %     type  >> help chemical_dataset, for more information
        err_type = 'MSE';
        wash_out = 0;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        [input_sequence,output_sequence] = chemical_dataset;
        t =  randperm(length(input_sequence),length(input_sequence));
        
        input_sequence = input_sequence(:,t)';
        output_sequence = output_sequence(:,t)';
        
        %% Regression
    case 'boston_housing'
        %     The modified Boston housing dataset consists of 489 data points, with each datapoint having 3 features. This dataset is a modified version of the Boston Housing dataset found on the UCI Machine Learning Repository.
        err_type = 'RMSE';
        wash_out = 0;
        config.train_fraction=0.6;    config.val_fraction=0.2;    config.test_fraction=0.2;
        
        data = load('boston.txt');
        
        input_sequence = data(:,1:3);
        output_sequence = data(:,4);
        
        t =  randperm(length(input_sequence),length(input_sequence));
        
        input_sequence = input_sequence(t,:);
        output_sequence = output_sequence(t,:);
        % rescale output
        [output_sequence] = featureNormailse(output_sequence,config);
        
        %% For specific methods and reservoirs
    case 'CPPN'
        err_type = 'empty';
        config.train_fraction=0.1;    config.val_fraction=0.1;    config.test_fraction=0.1;
        wash_out = 1;
        input_sequence= zeros(100,config.CPPN_inputs);
        output_sequence= zeros(100,config.CPPN_outputs);
        
        
    case 'autoencoder'
        
        err_type = 'NMSE';
        config.train_fraction=0.7;    config.val_fraction=0.15;    config.test_fraction=0.15;
        wash_out = 0;
        noise_factor = 0;
        
        type = 'image';
        
        data_length = 2000;
        switch(type)
            case 'digit'
                % image
                [XTrain, YTrain, XTest, YTest] = load_mnist('./Support files/other/Datasets/MNIST');
                input_sequence = [XTrain; XTest];
                input_sequence = input_sequence(1:data_length,:,:);
                input_sequence = reshape(input_sequence,size(input_sequence,1),size(input_sequence,2).^2);
                output_sequence = input_sequence;
                
                % imshow(reshape(input_sequence(8,:),28,28))
                
            case 'denoising_digit'
                % image
                [XTrain, YTrain, XTest, YTest] = load_mnist('./Support files/other/Datasets/MNIST');
                input_sequence = [XTrain; XTest];
                input_sequence = input_sequence(1:data_length,:,:);
                input_sequence = reshape(input_sequence,size(input_sequence,1),size(input_sequence,2).^2);
                output_sequence = input_sequence;
                % add noise
                input_sequence = input_sequence+ randn(size(output_sequence,1),1)*noise_factor;
                
            case 'image'
                
                [XTrain, ~,XTest] = load_cifar10;
                
                input_sequence = [XTrain; XTest];
                input_sequence = input_sequence(1:data_length,:,:);
                output_sequence = input_sequence;
                
                %                 I = imread('cat_50pixel.jpg');
                %                 output_sequence =double(I(:))./255;
                %                 input_sequence = double(I(:))./255 + randn(size(output_sequence))*noise_factor;
                %
        end
        
        %input_sequence= input_sequence';
        %output_sequence= output_sequence';
        
    case 'image_painting'
        
        err_type = 'NMSE';
        config.train_fraction=1;    config.val_fraction=0;    config.test_fraction=0;
        wash_out = 0;
        data_length = 1;
        
        %get image data
        I = imread('cat.jpg');
        
        %I = imread('panda400px.jpg');
        
        %         [XTrain, ~,XTest] = load_cifar10;
        %         I = [XTrain; XTest];
        %I = reshape(I(1:data_length,:),data_length,size(I,2)/3,3);
        
        % given input x, y coordinates
        [x, y] = meshgrid(1:1:size(I,1));
        input_sequence = [x(:),y(:)];
        %         [x, y] = meshgrid(1:1:sqrt(size(I,2)/3));
        %         input_sequence = repmat([x(:),y(:)],data_length,1);
        
        % recreate image
        r = double(I(:,:,1))./255;
        g = double(I(:,:,2))./255;
        b = double(I(:,:,3))./255;
        
        %         r = double(I(1:data_length,1:1024))./255;
        %         g = double(I(1:data_length,1025:2048))./255;
        %         b = double(I(1:data_length,2049:end))./255;
        
        output_sequence = [r(:) g(:) b(:)];
        
        %% Metrics
        case 'MC'
        
        err_type = 'MC';
        config.train_fraction=0.5;    config.val_fraction=0;    config.test_fraction=0.5;
        wash_out = 50;
        
        data_length = 1000 + wash_out*2;
        
        n_internal_units = config.total_units;%sum(config.num_nodes);
        
        n_output_units = n_internal_units*2;
        n_input_units = 1;
        
        data_sequence = 2*rand(n_input_units,data_length+1+n_output_units)-1;
        
        % rescale for each reservoir
        %[data_sequence] = featureNormailse(data_sequence,config);
        
        if config.discrete %strcmp(config.res_type,'elementary_CA') || strcmp(config.res_type,'2d_CA') || strcmp(config.res_type,'RBN')
            data_sequence = floor(heaviside(data_sequence));
        end
        
        input_sequence = data_sequence(n_output_units+1:data_length+n_output_units)';
        
        for i = 1:n_output_units
            output_sequence(:,i) = data_sequence(n_output_units+1-i:data_length+n_output_units-i);
        end              
end

if classification_data
    
    wash_out = 0;
    
    % make sure distribution is even
    indx = [];
    num_classes = size(output_sequence,2);
    div_d = ceil(size(input_sequence,1)/num_classes);
    for k = 1:num_classes
        n = find(output_sequence(:,k));
        indx(:,k) = n(randi([1 length(n)],1,div_d));
    end
    
    % get correct lengths and id's for data split
    train_id = indx(1:div_d*config.train_fraction,:);
    val_id = indx(1:div_d*config.val_fraction,:);
    test_id = indx(1:div_d*config.test_fraction,:);
    
    % randomise
    train_id = train_id(randperm(length(train_id(:))));
    val_id = val_id(randperm(length(val_id(:))));
    test_id = test_id(randperm(length(test_id(:))));
    
    input_sequence = input_sequence([train_id val_id test_id],:);
    output_sequence = output_sequence([train_id val_id test_id],:);
end

%% preprocessing
switch(config.input_mechanism)
    case'spiking'
        config.preprocess = 'scaling'; % must be between [0 1]
        config.preprocess_shift = [0 1];
        % rescale training data
        [input_sequence] = featureNormailse(input_sequence,config);
        % find appropriate filter for data
        config.num_gens =2000;
        config.max_order=48;
        config.max_period = 3;
        
        [config.filter] = filterGA(input_sequence,config);
        wash_out =0;
    otherwise
        % rescale training data
        [input_sequence] = featureNormailse(input_sequence,config);
        
end

% split datasets
[train_input_sequence,val_input_sequence,test_input_sequence] = ...
    split_train_test3way(input_sequence,config.train_fraction,config.val_fraction,config.test_fraction);

[train_output_sequence,val_output_sequence,test_output_sequence] = ...
    split_train_test3way(output_sequence,config.train_fraction,config.val_fraction,config.test_fraction);

% Add extra for washout
%train_input_sequence= [train_input_sequence(1:wash_out,:); train_input_sequence];
%train_output_sequence= [train_output_sequence(1:wash_out,:); train_output_sequence];

if config.val_fraction > 0
    if (size(val_input_sequence,1) < wash_out)
        wash_input_sequence = [];
        wash_output_sequence = [];
        while(size(wash_input_sequence,1) < size(val_input_sequence,1) + wash_out)
            wash_input_sequence = [val_input_sequence; wash_input_sequence];
            wash_output_sequence = [val_output_sequence; wash_output_sequence];
        end
        val_input_sequence = wash_input_sequence;
        val_output_sequence = wash_output_sequence;
    else
        %val_input_sequence= [val_input_sequence(1:wash_out,:); val_input_sequence];
        %val_output_sequence= [val_output_sequence(1:wash_out,:); val_output_sequence];
    end
end

if config.test_fraction > 0
    if size(test_input_sequence,1)< wash_out
        wash_input_sequence = [];
        wash_output_sequence = [];
        while(size(wash_input_sequence,1) < size(test_input_sequence,1) + wash_out)
            wash_input_sequence = [test_input_sequence; wash_input_sequence];
            wash_output_sequence = [test_output_sequence; wash_output_sequence];
        end
        test_input_sequence = wash_input_sequence;
        test_output_sequence = wash_output_sequence;
    else
        %test_input_sequence= [test_input_sequence(1:wash_out,:); test_input_sequence];
        %test_output_sequence= [test_output_sequence(1:wash_out,:); test_output_sequence];
    end
end

% squash into structure
config.train_input_sequence = train_input_sequence;
config.train_output_sequence = train_output_sequence;
config.val_input_sequence = val_input_sequence;
config.val_output_sequence = val_output_sequence;
config.test_input_sequence = test_input_sequence;
config.test_output_sequence = test_output_sequence;

config.wash_out = wash_out;
config.err_type = err_type;

% % if multi-objective, update input/output units
if ~isfield(config,'nsga2')
    config.task_num_inputs = size(config.train_input_sequence,2);
    config.task_num_outputs = size(config.train_output_sequence,2);
end

% Go back to old seed
rng(temp_seed,'twister');
