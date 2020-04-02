


%%Time specifications:
Fs = 8000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.1;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc = 30;                     % hertz
data = (sin(2*pi*Fc*t)+1)/2;

%data = rand(1000,1);

config.dataset = 'laser';          % Task to evolve for
config.num_nodes = 1;
config.res_type = '';

% get any additional params. This might include:
% details on reservoir structure, extra task variables, etc.
config = getAdditionalParameters(config);
config.preprocess = 'scaling';
config.prepocess_shift = 'zero to one';
        
% get dataset information
config = selectDataset(config);
data = config.train_input_sequence;


%% run GA
config.num_gens =1000;
config.max_order=48;
config.max_period = 2;
config.mutate_type = 'gaussian';

[filter] = filterGA(data(:,1),config);
