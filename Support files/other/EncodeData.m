function varargout = EncodeData(varargin)

data_sequence = varargin{1};
config = varargin{2};
if nargin >= 3
    output_sequence = varargin{3};
end

%% preprocessing
switch(config.input_mechanism)
    case 'BSA'
        config.preprocess = 'scaling'; % must be between [0 1]
        config.preprocess_shift = [0 1];
        % rescale training data
        [data_sequence] = featureNormailse(data_sequence,config);
        % find appropriate filter for data
        config.num_gens =2000;
        config.max_order=48;
        config.max_period = 3;
        
        [config.filter] = filterGA(data_sequence,config);
        %config.wash_out =0;
        
    case 'spike-rate'
        
        % define spike periods
        %         config.Tmin = 1;
        %         config.Tmax = 20;
        N = length(data_sequence);
        
        % normalise input signal s(t) to s_hat(t)
        config.preprocess = 'scaling';
        config.preprocess_shift = [0 1]; % range for data
        [norm_input_sequence] = featureNormailse(data_sequence,config);
        
        % plot
        subplot(2,1,1)
        plot(norm_input_sequence)
        
        Sv = zeros(N,1);
        
        i0 = 1;
        i1 = 1;
        
        while (i0 < N)
            spike_period = floor(norm_input_sequence(i0)*(config.Tmax-config.Tmin)) + config.Tmin;
            i1 = i0+spike_period;
            Sv(i1) = 1;
            i0 = i1;
        end
        
        subplot(2,1,2)
        plot(Sv)
        
        data_sequence = Sv(1:N);
        
    otherwise
        
        % rescale training data
        [data_sequence] = featureNormailse(data_sequence,config);
        if nargin == 3
            [output_sequence] = featureNormailse(output_sequence,config);
        end
        
end

varargout{1} = data_sequence;
varargout{2} = config;
if nargin == 3
    varargout{3} = output_sequence;
else
    varargout{3} = [];
end