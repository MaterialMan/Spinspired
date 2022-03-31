% Separation Metrics and Kernel Quality
function [metrics,M] = getMetrics(individual,config)

scurr = rng;
temp_seed = scurr.Seed;

% set parameters
metrics = [];
config.reg_param = 10e-9;
config.wash_out = 50;
metrics_type =  config.metrics;
num_timesteps = round(individual.total_units*1.5) + config.wash_out; % input should be twice the size of network + wash out
MC_num_timesteps = 500 + config.wash_out*2;
n_input_units = individual.n_input_units;

for metric_item = 1:length(config.metrics)
    
    rng(1,'twister');
    clearvars -except metrics num_timesteps MC_num_timesteps n_input_units individual config metric_item metrics_type temp_seed
    
    switch metrics_type{metric_item}
        
        % kernel rank
        case 'KR'
            
            %define input signal
            ui = 2*rand(num_timesteps,n_input_units)-1;
            
            input_sequence = repmat(ui(:,1),1,n_input_units);
            
            % rescale for each reservoir
            %[input_sequence] = featureNormailse(input_sequence,config);
            
            %kernel matrix - pick 'to' at halfway point
            M = config.assessFcn(individual,input_sequence,config);
            
            %catch errors
            M(isnan(M)) = 0;
            M(isinf(M)) = 0;
            
            %% Kernal Quality
            s = svd(M);
            
            tmp_rank_sum = 0;
            full_rank_sum = 0;
            e_rank = 1;
            for i = 1:length(s)
                full_rank_sum = full_rank_sum + s(i);
                while (tmp_rank_sum < full_rank_sum * 0.99)
                    tmp_rank_sum = tmp_rank_sum + s(e_rank);
                    e_rank= e_rank+1;
                end
            end
            
            kernel_rank = e_rank-1;
                   
            metrics = [metrics kernel_rank];
            
            % Genralization Rank
        case 'GR'
            % define input signal
            input_sequence = 0 + 0.1*rand(num_timesteps,n_input_units)-0.05;
            
            % rescale for each reservoir
            %[input_sequence] = featureNormailse(input_sequence,config);
            
            %collect states
            G = config.assessFcn(individual,input_sequence,config);
            
            %catch errors
            G(isnan(G)) = 0;
            G(isinf(G)) = 0;
            
            % get rank of matrix
            s = svd(G);
            
            %claculate effective rank
            tmp_rank_sum = 0;
            full_rank_sum = 0;
            e_rank = 1;
            for i = 1:length(s)
                full_rank_sum = full_rank_sum +s(i);
                while (tmp_rank_sum < full_rank_sum * 0.99)
                    tmp_rank_sum = tmp_rank_sum + s(e_rank);
                    e_rank= e_rank+1;
                end
            end
            gen_rank = e_rank-1;
            
            metrics = [metrics gen_rank];
            
            % LE measure
        case 'LE'
            seed = 1;
            LE = lyapunovExponent(individual,config,seed);
            metrics = [metrics LE];
            
            % Entropy measure
        case 'entropy'
            
            data_length = num_timesteps;%individual.total_units*2 + config.wash_out;%400;
            input_sequence = ones(data_length,n_input_units);
            [input_sequence] = featureNormailse(input_sequence,config);
            
            X = config.assessFcn(individual,input_sequence,config);
            C = X'*X;
            
            X_eig = eig(C);
            
            normX_eig = X_eig./sum(X_eig);
            
            H = -sum(normX_eig.*log2(normX_eig));
            
            entropy = real(H/log2(size(X,2)));
            
            entropy(isnan(entropy)) = 0;
            metrics = [metrics entropy*100];
            
            % linear memory capacity
        case 'linearMC'
            
            % measure MC multiple times
            mc_seed = 1;
            temp_MC = testMC(individual,config,mc_seed,MC_num_timesteps);
            
            MC = mean(temp_MC);
            
            metrics = [metrics MC];
            
            % quadratic memory capacity (nonlinear) 
        case 'quadMC'
            
            quad_MC = quadraticMC(individual,config,1);
            
            metrics = [metrics quad_MC];
            
            % cross-memory capacity (nonlinear) 
        case 'crossMC'
            
            cross_MC = crossMC(individual,config,1);
            
            metrics = [metrics cross_MC];
            
            % separation property
        case 'separation'
            
            data_length = num_timesteps;%individual.total_units*4 + config.wash_out*2;%400;
            
            u1 = (rand(data_length,n_input_units)-1);
            u2 = (rand(data_length,n_input_units));
            [u1] = featureNormailse(u1,config);
            [u2] = featureNormailse(u2,config);
            
            D= norm(u1-u2);
            
            X1 = config.assessFcn(individual,u1,config);
            
            X2 = config.assessFcn(individual,u2,config);
            
            sep = norm(X1 - X2)/D;
            
            metrics = [metrics sep];
            %abs(X1-X2)/D(i,config.wash_out+1:end);
            
            
            %             input_sequence = ones(data_length,1);
            %
            %             X = config.assessFcn(individual,input_sequence,config);
            %
            %             centre_of_mass = mean(X);
            %
            %             inter_class_distance =
            %
            %             intra_class_var =
            %
            %             sep = inter_class_distance/(intra_class_var + 1);
            %
            
        case 'mutalInformation'
            
            data_length = individual.total_units*4 + config.wash_out*2;%400;
            
            input_sequence = (rand(data_length,n_input_units)-1).*config.scaler;
            
            [input_sequence] = featureNormailse(input_sequence,config);
            
            X = config.assessFcn(individual,input_sequence,config);
            
            for i = 1:size(X,1)
                for j = 1:size(X,1)
                    MI(j) = mutualInformation(X(i+1,:), X(i,j));
                end
                meanMI = mean(MI);
            end
            
        case 'transferEntropy'
            TE = transferEntropy(X, Y, W, varargin);
            
        case 'connectivity'
            metrics = [metrics individual.connectivity*(individual.total_units+ (config.add_input_states*config.task_num_inputs))];
            
        case 'matrix_dissimilarity'
            
            W = []; W_i =[];
            for i = 1:size(individual.W,1)
                for j = 1:size(individual.W,2)
                    W_i = [W_i individual.W{i,j}.*individual.W_scaling(i,j)];
                end
                W = [W; W_i];
            end
            
            dissimilarity = (norm(zeros(size(individual.W,1))-W)./20)*(individual.total_units + (config.add_input_states*config.task_num_inputs));
            
            metrics = [metrics dissimilarity];
    end
end

rng(temp_seed,'twister');