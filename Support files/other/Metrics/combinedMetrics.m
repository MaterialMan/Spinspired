function metrics = combinedMetrics(individual,config,seed,data_length)

rng(seed,'twister');

metrics = [];

%% Run network with random input
% define network sizes
n_internal_units = individual.total_units;
n_output_units = n_internal_units*2;
n_input_units = individual.n_input_units;

% define washout period
config.wash_out = n_output_units;

%define input signal
data_sequence = 2*rand(data_length,n_input_units)-1;
if config.discrete %strcmp(config.res_type,'elementary_CA') || strcmp(config.res_type,'2d_CA') || strcmp(config.res_type,'RBN')
    data_sequence = floor(heaviside(data_sequence));
end
input_sequence = repmat(data_sequence(:,1),1,n_input_units);

%kernel matrix - pick 'to' at halfway point
states = config.assessFcn(individual,input_sequence,config);

%catch errors
states(isnan(states)) = 0;
states(isinf(states)) = 0;

% split states
sequence_length = data_length/2;
train_states = states(1:sequence_length,:);
test_states = states(1+sequence_length:end,:);

% Find best reg parameter
reg_param = [10e-7];
config.err_type = 'NMSE';

%% define input data
for n_metrics = 1:length(config.metrics_names)
    
    switch(config.metrics_names{n_metrics})
        
        case 'L-MC'
            try load('L_MCdata.mat'); 
            catch mem_input_sequence = zeros(size(input_sequence(config.wash_out+1:data_length)));
            end
            if ~isequal(mem_input_sequence,input_sequence(config.wash_out+1:data_length))
               
                mem_input_sequence = input_sequence(config.wash_out+1:data_length);
                
                for i = 1:n_output_units
                    mem_output_sequence(:,i) = input_sequence(config.wash_out+1-i:data_length-i);
                end
                save('L_MCdata.mat','mem_input_sequence','mem_output_sequence');
            end
            
        case 'Q-MC'
            try load('Q_MCdata.mat'); 
            catch quad_input_sequence = zeros(size(input_sequence(config.wash_out+1:data_length)));
            end
            if ~isequal(quad_input_sequence,input_sequence(config.wash_out+1:data_length))
                quad_input_sequence = input_sequence(config.wash_out+1:data_length);
                for k = 1:n_output_units
                    for n = 1:length(quad_input_sequence)
                        if n > k
                            quad_output_sequence(n,k) = (3.*(quad_input_sequence(n-k,1).^2)) - 1;
                        else
                            quad_output_sequence(n,k) = quad_input_sequence(n,1);
                        end
                    end
                end
                save('Q_MCdata.mat','quad_input_sequence','quad_output_sequence');
            end
        case 'C-MC'
            try load('C_MCdata.mat'); 
            catch cross_input_sequence = zeros(size(input_sequence(config.wash_out+1:data_length)));
            end
            if ~isequal(cross_input_sequence,input_sequence(config.wash_out+1:data_length))
                cross_input_sequence = input_sequence(config.wash_out+1:data_length);
                for n = 1:length(cross_input_sequence)
                    cnt = 1;
                    for k = 1:n_output_units
                        for k_p = 1:n_output_units
                            if (n > k && n > k_p) && k ~= k_p
                                %crossFcn(k,k_p) = cross_input_sequence(n-k,1)*cross_input_sequence(n-k_p,1);
                                cross_output_sequence(n,cnt) = cross_input_sequence(n-k,1)*cross_input_sequence(n-k_p,1);
                            else
                                %crossFcn(k,k_p) = 0;
                                cross_output_sequence(n,cnt) = 0;
                            end
                            cnt = cnt+1;
                        end
                        %cross_output_sequence(n,k) = sum(crossFcn(k,:));
                    end
                    % cross_output_sequence(n,:) = crossFcn(:);
                end
                save('C_MCdata.mat','cross_input_sequence','cross_output_sequence');
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
            
        case 'L-MC'
            %% calculate linear memory capacity (MC)            
            train_input_sequence = mem_input_sequence(1:sequence_length,:);
            test_input_sequence = mem_input_sequence(1+sequence_length:end,:);
            
            train_output_sequence = mem_output_sequence(1:sequence_length,:);
            test_output_sequence = mem_output_sequence(1+sequence_length:end,:);
            
            % get output and weights
            Y = calculateWeights(train_states,train_output_sequence,test_states,test_output_sequence,reg_param,config);
            
            MC_k= 0; Cm = 0;
            test_in_var = test_input_sequence(:,1);
            targVar = 1/(length(test_in_var)-1) * sum((test_in_var-mean(test_in_var)).*(test_in_var-mean(test_in_var)));
            
            for i = 1:n_output_units
                
                coVar = 1/(length(Y(:,i))-1) * sum((test_output_sequence(:,i)-mean(test_output_sequence(:,i)))...
                    .*(Y(:,i)-mean(Y(:,i))));
                outVar = 1/(length(Y(:,i))-1) * sum((Y(:,i)-mean(Y(:,i))).*(Y(:,i)-mean(Y(:,i))));
                totVar = (outVar*targVar(1));
                MC_k(i) = (coVar*coVar)/totVar;
                
                %remove low values from measure
                if MC_k(i) <= 0.1
                    MC_k(i) = 0;
                end
            end
            
            MC_k(isnan(MC_k)) = 0;
            MC = sum(MC_k);
            
            metrics = [metrics MC];
            
        case 'Q-MC'
            %% %% quadratic memory capacity - y_k(n)=3u^2(n-k)-1       
            train_input_sequence = quad_input_sequence(1:sequence_length,:);
            test_input_sequence = quad_input_sequence(1+sequence_length:end,:);
            
            train_output_sequence = quad_output_sequence(1:sequence_length,:);
            test_output_sequence = quad_output_sequence(1+sequence_length:end,:);
            
            % get output and weights
            Y = calculateWeights(train_states,train_output_sequence,test_states,test_output_sequence,reg_param,config);
            
            % calculate output
            for i = 1:size(Y,2)
                
                nmse = (mean((test_output_sequence(:,i)-Y(:,i)).^2)/mean(test_output_sequence(:,i).^2));
                
                if nmse > 1
                    nmse = 1;
                end
                
                quad_MC(i) = 1 - nmse;
                
                %remove low values from measure
                if quad_MC(i) <= 0.1
                    quad_MC(i) = 0;
                end
                
            end
            quad_MC = sum(quad_MC);
            
            metrics = [metrics quad_MC];
            
        case 'C-MC'
            %% cross memory capacity - y_kk'(n)=u(n-k)u(n-k')
            train_input_sequence = cross_input_sequence(1:sequence_length,:);
            test_input_sequence = cross_input_sequence(1+sequence_length:end,:);
            
            train_output_sequence = cross_output_sequence(1:sequence_length,:);
            test_output_sequence = cross_output_sequence(1+sequence_length:end,:);
            
            % get output and weights
            Y = calculateWeights(train_states,train_output_sequence,test_states,test_output_sequence,reg_param,config);
            
            % calculate output
            for i = 1:size(Y,2)
                
                if sum(test_output_sequence(:,i)) > 0
                    nmse = (mean((test_output_sequence(:,i)-Y(:,i)).^2)/mean(test_output_sequence(:,i).^2));
                    
                    if nmse > 1
                        nmse = 1;
                    end
                    
                    cross_MC(i) = 1 - nmse;
                    
                    %remove low values from measure
                    if cross_MC(i) <= 0.1
                        cross_MC(i) = 0;
                    end
                    
                else
                    cross_MC(i) = 0;
                end
                
            end
            cross_MC = sum(cross_MC);
            
            metrics = [metrics cross_MC];
            
        otherwise
            
    end
end
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