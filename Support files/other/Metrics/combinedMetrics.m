%% function assesses all given metrics using a single data input 

function metrics = combinedMetrics(individual,config,seed,data_length)

rng(seed,'twister');

metrics = [];
config.err_type = 'NRMSE';

%% Run network with random input
% define network sizes
n_internal_units = individual.total_units;
n_output_units = n_internal_units*2;
n_input_units = individual.n_input_units;

% define washout period
config.wash_out = n_output_units;

%define input signal
data_sequence = 2*rand(data_length+config.wash_out,n_input_units)-1;
if config.discrete %strcmp(config.res_type,'elementary_CA') || strcmp(config.res_type,'2d_CA') || strcmp(config.res_type,'RBN')
    data_sequence = floor(heaviside(data_sequence));
end
input_sequence = repmat(data_sequence(:,1),1,n_input_units);

% no preprocessing, i.e. rescaling
config.preprocess = ''; 
config.preprocess_shift = [];
[input_sequence, config] = EncodeData(input_sequence,config);
            
%kernel matrix - pick 'to' at halfway point
states = config.assessFcn(individual,input_sequence,config);

%catch errors
states(isnan(states)) = 0;
states(isinf(states)) = 0;

% split states
sequence_length = data_length/2;
data_length = length(data_sequence);
train_states = states(1:sequence_length,:);
test_states = states(1+sequence_length:end,:);

% Find best reg parameter
reg_param = [10e-7];
%config.err_type = 'NMSE';

%% define input data
for n_metrics = 1:length(config.metrics_names)
    
    switch(config.metrics_names{n_metrics})
        
        case 'L_MC'
            try load(strcat('L_MCdata',num2str(individual.total_units),'.mat')); 
            catch mem_input_sequence = zeros(size(input_sequence(config.wash_out+1:data_length)));
            end
            if ~isequal(mem_input_sequence,input_sequence(config.wash_out+1:data_length))
               
                mem_input_sequence = input_sequence(config.wash_out+1:data_length);
                
                for i = 1:n_internal_units%n_output_units
                    mem_output_sequence(:,i) = input_sequence(config.wash_out+1-i:data_length-i);
                end
                save(strcat('L_MCdata',num2str(individual.total_units),'.mat'),'mem_input_sequence','mem_output_sequence');
            end
            
        case 'Q_MC' %Quadratic memory: yk(n)=(3u^2(n-k)-1)/2  ----- taken from Duport et al. and Ortin and Pesquera (2019, doi.org/10.3389/fphy.2019.00210)
            try load(strcat('Q_MCdata',num2str(individual.total_units),'.mat')); 
            catch quad_input_sequence = zeros(size(input_sequence(config.wash_out+1:data_length)));
            end
            if ~isequal(quad_input_sequence,input_sequence(config.wash_out+1:data_length))
                quad_input_sequence = input_sequence(config.wash_out+1:data_length);
                for i = 1:n_internal_units%n_output_units%n_internal_units
                    quad_output_sequence(:,i) = ((input_sequence(config.wash_out+1-i:data_length-i).^2)*3 -1)/2;
                end
                
                save(strcat('Q_MCdata',num2str(individual.total_units),'.mat'),'quad_input_sequence','quad_output_sequence');
            end
            
        case 'C_MC' %Cubic memory: yk(n)=(5u^3(n-k)-3u(n-k))/2
            try load(strcat('C_MCdata',num2str(individual.total_units),'.mat'));
            catch cubic_input_sequence = zeros(size(input_sequence(config.wash_out+1:data_length)));
            end
            if ~isequal(cubic_input_sequence,input_sequence(config.wash_out+1:data_length))
                cubic_input_sequence = input_sequence(config.wash_out+1:data_length);
                for i = 1:n_internal_units%n_output_units
                    cubic_output_sequence(:,i) =  (5 * (input_sequence(config.wash_out+1-i:data_length-i).^3) - (3 * input_sequence(config.wash_out+1-i:data_length-i)))/2;
                end
                
                save(strcat('C_MCdata',num2str(individual.total_units),'.mat'),'cubic_input_sequence','cubic_output_sequence');
            end
            
        case 'X_MC' %Cross memory: yk,k'=u(n-k)·u(n-k') summing over k,k' for k < k'
            try load(strcat('X_MCdata',num2str(individual.total_units),'.mat')); 
            catch cross_input_sequence = zeros(size(input_sequence(config.wash_out+1:data_length)));
            end
            if ~isequal(cross_input_sequence,input_sequence(config.wash_out+1:data_length))
                cross_input_sequence = input_sequence(config.wash_out+1:data_length);
                cnt = 1;
                for k = 1:n_internal_units
                    for k_p = 1:n_internal_units%n_output_units
                        if k < k_p
                            cross_output_sequence(:,cnt) = input_sequence(config.wash_out+1-k:data_length-k,1).*input_sequence(config.wash_out+1-k_p:data_length-k_p,1);
                            cnt = cnt+1;
                        end
                    end
                end
                save(strcat('X_MCdata',num2str(individual.total_units),'.mat'),'cross_input_sequence','cross_output_sequence');
            end
        otherwise
    end
end

%% assess all metrics
for n_metrics = 1:length(config.metrics_names)
    
    switch(config.metrics_names{n_metrics})
        
        case 'KR'
            %% Kernal Quality
            s = svd(states);
            
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
            
        case 'L_MC'
            %% calculate linear memory capacity (MC)                       
            C_y = calculateCapacity(mem_output_sequence(1:sequence_length,:),mem_output_sequence(1+sequence_length:end,:)...
                ,train_states,test_states,reg_param,config);
            
            metrics = [metrics C_y];            
            
        case 'Q_MC'
            %% %% quadratic memory capacity - (y_k(n)= 3u^2(n-k)-1)/2       
            C_y = calculateCapacity(quad_output_sequence(1:sequence_length,:),quad_output_sequence(1+sequence_length:end,:)...
                ,train_states,test_states,reg_param,config);
            
            metrics = [metrics C_y];
            
        case 'C_MC'
            %% %% cubic memory capacity - yk(n)=(5u^3(n-k)-3u(n-k))/2     
            C_y = calculateCapacity(cubic_output_sequence(1:sequence_length,:),cubic_output_sequence(1+sequence_length:end,:)...
                ,train_states,test_states,reg_param,config);
            
            metrics = [metrics C_y];
            
        case 'X_MC'
            %% cross memory capacity - y_kk'(n)=u(n-k)u(n-k')
            C_y = calculateCapacity(cross_output_sequence(1:sequence_length,:),cross_output_sequence(1+sequence_length:end,:)...
                ,train_states,test_states,reg_param,config);
            
            metrics = [metrics C_y];
            
        otherwise
            
    end
end
end

function C_y = calculateCapacity(train_output_sequence,test_output_sequence,train_states,test_states,reg_param,config)

Y = calculateWeights(train_states,train_output_sequence,test_states,test_output_sequence,reg_param,config);

% calculate output
for i = 1:size(Y,2)
    
    C_yk(i) = 1 - (sum((Y(:,i)-test_output_sequence(:,i)).^2)/sum(test_output_sequence(:,i).^2));
    
    %remove low values from measure
    if C_yk(i) < 0.1
        C_yk(i) = 0;
    end
    
end

C_y = sum(C_yk);

end

function [Y,output_train_sequence,output_test_sequence,output_weights] = calculateWeights(train_states,train_output_sequence,test_states,test_output_sequence,reg_param,config)

reg_train_error = [];
reg_weights=[];

for i = 1:length(reg_param)
    %Train: tanspose is inversed compared to equation
    output_weights = train_output_sequence'*train_states*inv(train_states'*train_states + reg_param(i)*eye(size(train_states'*train_states)));
    
    % Calculate trained output Y
    output_train_sequence = train_states*output_weights';
    reg_train_error(i,:)  = calculateError(output_train_sequence,train_output_sequence,config);
    
    % Calculate trained output Y
    %output_test_sequence = test_states*output_weights';
    %reg_test_error(i,:)  = calculateError(output_test_sequence,test_output_sequence,config);
    reg_weights(i,:,:) =output_weights;
end

[~, reg_indx]= min(sum(reg_train_error,2));
output_weights =reshape(reg_weights(reg_indx,:,:),size(reg_weights,2),size(reg_weights,3));

%remove NaNs
output_weights(isnan(output_weights)) = 0;

if config.discrete
    Y = double((test_states * output_weights') > 0);
else
    Y = test_states * output_weights';
end

end