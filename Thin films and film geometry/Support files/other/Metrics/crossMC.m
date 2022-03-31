function cross_MC = crossMC(individual,config,seed)

rng(seed,'twister');

n_internal_units = individual.total_units;%sum([genotype.nInternalUnits]);

n_output_units = n_internal_units*2;
n_input_units = individual.n_input_units;

%% Assign input data and collect target output
data_length = n_output_units*3 + config.wash_out*2;
sequence_length = (data_length-n_output_units)/2;%200;
data_sequence = 2*rand(n_input_units,data_length)-1;

% rescale for each reservoir
%[data_sequence] = featureNormailse(data_sequence,config);

if config.discrete %strcmp(config.res_type,'elementary_CA') || strcmp(config.res_type,'2d_CA') || strcmp(config.res_type,'RBN')
    data_sequence = floor(heaviside(data_sequence));
end

input_sequence = data_sequence';
%input_sequence = data_sequence(n_output_units+1:data_length+n_output_units)';

%% cross memory capacity - y_kk'(n)=u(n-k)u(n-k')
cross_input_sequence = input_sequence;
for n = 1:length(input_sequence)
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

% config.preprocess = 'rescale';
% [cross_output_sequence] = featureNormailse(cross_output_sequence,config);
cross_input_sequence = cross_input_sequence(n_output_units+1:end,:);
cross_output_sequence = cross_output_sequence(n_output_units+1:end,:);

%cross_output_sequence(cross_output_sequence == 0,:) = [];

% split data
cross_train_input_sequence = repmat(cross_input_sequence(1:sequence_length,:),1,n_input_units);
cross_test_input_sequence = repmat(cross_input_sequence(1+sequence_length:end,:),1,n_input_units);

cross_train_output_sequence = cross_output_sequence(1:sequence_length,:);
cross_test_output_sequence = cross_output_sequence(1+sequence_length:end,:);

%train
cross_states_train = config.assessFcn(individual,cross_train_input_sequence,config);
cross_output_weights = cross_train_output_sequence(config.wash_out+1:end,:)'*cross_states_train*inv(cross_states_train'*cross_states_train + config.reg_param*eye(size(cross_states_train'*cross_states_train)));
cross_Y_train = cross_states_train * cross_output_weights';

% test
cross_states_test =  config.assessFcn(individual,cross_test_input_sequence,config);
cross_Y_test =cross_states_test * cross_output_weights';

% calculate output
for i = 1:size(cross_Y_test,2)
    
    %nmse = (mean((cross_test_output_sequence(config.wash_out+1:end,i)-cross_Y_test(:,i)).^2)/var(cross_test_output_sequence(config.wash_out+1:end,i)));
    %nmse = (mean((cross_train_output_sequence(config.wash_out+1:end,i)-cross_Y_train(:,i)).^2)/var(cross_train_output_sequence(config.wash_out+1:end,i)));
    if sum(cross_test_output_sequence(config.wash_out+1:end,i)) > 0
        nmse = (mean((cross_test_output_sequence(config.wash_out+1:end,i)-cross_Y_test(:,i)).^2)/mean(cross_test_output_sequence(config.wash_out+1:end,i).^2));
        
        %    if sum(cross_Y_test(:,i)) == 0 || isnan(sum(cross_Y_test(:,i))) || isinf(sum(cross_Y_test(:,i)))
        %        nmse = 1;
        %    else
        %     mse=norm(cross_test_output_sequence(config.wash_out+1:end,i)-cross_Y_test(:,i),2)^2/length(cross_test_output_sequence(config.wash_out+1:end,i));
        %     sigEner=norm(cross_test_output_sequence(config.wash_out+1:end,i))^2;
        %     nmse=(mse/sigEner);
        %
        
        %    end
        %
        
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

% cross_MC_k= 0;
% test_in_var = cross_test_input_sequence(config.wash_out+1:end,1);
% targVar = 1/(length(test_in_var)-1) * sum((test_in_var-mean(test_in_var)).*(test_in_var-mean(test_in_var)));
%
% for i = 1:n_output_units.^2
%
%     coVar = 1/(length(cross_Y(:,i))-1) * sum((cross_test_output_sequence(config.wash_out+1:end,i)-mean(cross_test_output_sequence(config.wash_out+1:end,i)))...
%        .*(cross_Y(:,i)-mean(cross_Y(:,i))));
%     outVar = 1/(length(cross_Y(:,i))-1) * sum((cross_Y(:,i)-mean(cross_Y(:,i))).*(cross_Y(:,i)-mean(cross_Y(:,i))));
%     totVar = (outVar*targVar(1));
%     cross_MC_k(i) = (coVar*coVar)/totVar;
%
% end
%cross_MC_k(isnan(cross_MC_k)) = 0;
%cross_MC = sum(cross_MC_k);

end