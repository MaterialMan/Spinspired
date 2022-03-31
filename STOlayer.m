classdef STOlayer < nnet.layer.Layer & nnet.layer.Formattable

    properties
        % Layer properties.
        NumHiddenUnits
        OutputMode
        %HiddenState
        % extras
        config
        population
    end

    properties (Learnable)
        % Layer learnable parameters.
%         input_scaling
        input_weights
%         leak_rate
%         internal_scaling
 %        internal_weights
    end

%     properties (State)
%         % Layer state parameters.
%         HiddenState
%     end
    
    methods
        function layer = STOlayer(numHiddenUnits,inputSize,args)

            arguments
                numHiddenUnits
                inputSize
                args.Name = "";
                args.OutputMode = "sequence"
            end
    
            % type of network to evolve
            config.res_type = 'STO';            % state type of reservoir(s) to use. E.g. 'RoR' (Reservoir-of-reservoirs/ESNs), 'ELM' (Extreme learning machine), 'Graph' (graph network with multiple functions), 'DL' (delay line reservoir) etc. Check 'selectReservoirType.m' for more. Place reservoirs in cell ({}) for heterotic systems.
            config.num_nodes = {numHiddenUnits};                   % num of nodes in each sub-reservoir, e.g. if config.num_nodes = [10,5,15], there would be 3 sub-reservoirs with 10, 5 and 15 nodes each.
            config = selectReservoirType(config);         % collect function pointers for the selected reservoir type
            config.pop_size = 1;                       % initail population size. Note: this will generally bias the search to elitism (small) or diversity (large)
            config.dataset = '';
            config.test = 1;
            config.train_input_sequence = zeros(inputSize);
            config.train_output_sequence = zeros(inputSize);
            config.task_num_inputs = inputSize;
            config.task_num_outputs = inputSize;

            config = getAdditionalParameters(config);
            population = config.createFcn(config);

            % Set layer name.
            layer.NumHiddenUnits = numHiddenUnits;
            layer.Name = args.Name;
            layer.OutputMode = args.OutputMode;

            % Set layer description.
            layer.Description = "STO with " + numHiddenUnits + " channels";
        
            % washout
            %layer.wash_out = 50;

            % Initialize scaling params.
%             layer.internal_scaling = population.init_W_scaling;
%             layer.input_scaling = population.init_input_scaling;
%             layer.leak_rate = population.leak_rate;
% 
%             % input weights
             layer.input_weights = full(population.layer(1).input_weights{1}'); 
% 
%             % internal weights
%             layer.internal_weights = full(population.W);

            %layer.input_delay = population.input_delay;
            %layer.interpolation_length = population.interpolation_length;

            layer.config = config;
            layer.population = population;

            % Initialize layer states.
%            layer = resetState(layer);

        end

        function [Z] = predict(layer,X)
            % Forward input data through the layer at prediction time and
            % output the result and updated state.

            % Initialize sequence output.
            numHiddenUnits = layer.NumHiddenUnits;
            miniBatchSize = size(X,finddim(X,"B"));
            numTimeSteps = size(X,finddim(X,"T"));
            inputSize = size(X,finddim(X,"C"));

            if layer.OutputMode == "sequence"
                Z = zeros(numHiddenUnits+inputSize,miniBatchSize,numTimeSteps,"like",X);
                Z = dlarray(Z,"CBT");
            else
                Z = zeros(numHiddenUnits+inputSize,miniBatchSize,"like",X);
                Z = dlarray(Z,"CB");
            end

            % Calculate WX + b.
            X = stripdims(X);
            input_sequence = X;

            % reset params - i.e. params being learnt
            %layer.population.input_scaling = repmat(layer.input_scaling,1,numHiddenUnits);
            layer.population.layer(1).input_weights{1} = layer.input_weights';
            %layer.population.leak_rate = layer.leak_rate;
            %layer.population.W_scaling = double(extractdata(repmat(layer.internal_scaling,numHiddenUnits,numHiddenUnits)));
            %layer.population.W = double(extractdata(layer.internal_weights));

            % run system
            layer.config.wash_out = 0;

            for batch = 1:miniBatchSize

                if length(size(input_sequence)) > 2
                    input = reshape(input_sequence(:,batch,:),inputSize,numTimeSteps)'; 
                else
                    input = input_sequence';
                end

                %input_scaled = featureNormailse(input,layer.config);
                input_reform = double(extractdata(input));
                HiddenState = layer.config.assessFcn(layer.population,input_reform,layer.config);          
                HiddenState = [HiddenState input]; % must add input to keep dlarray tracing

                % store data
                if layer.OutputMode == "last"
                    Z(:,batch) = HiddenState(end,:);
                else
                    Z(:,batch,:) = HiddenState';
                end
            end
        end

%         function layer = resetState(layer)
%             % Reset layer state.
%             numHiddenUnits = layer.NumHiddenUnits;
%             layer.HiddenState = zeros(numHiddenUnits,1);
% 
%         end

    end
end