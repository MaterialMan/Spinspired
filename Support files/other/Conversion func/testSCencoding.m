%% Stocastic Computing encoding example
% converts continious number into a random bit stream, using the continious number as a probability 
% - larger N value will give greater accuarcy
original_input = rand(1,100);
N = 100;

%% encode
%  Preallocate the sc_array for performance
output = zeros(size(original_input, 1), N*size(original_input, 2));

for i = 1:size(original_input,1)
    for j = 1:size(original_input, 2)
        offset_j = (j-1)*N + 1;
        output(i,offset_j:offset_j+N-1) = rand(1,N)<original_input(i,j);
    end
end



%% fcn transformation
%fcn_output = tanh(output);

%% decode
% Determine the dimensions of the decimal matrix
converted_input = zeros(size(output, 1), size(output, 2) / N);

for i = 1:size(converted_input,1)
    for j = 1:size(converted_input, 2)
        offset_j = (j-1)*N + 1;
        converted_input(i,j) = sum(output(i,offset_j:offset_j+N-1)) / N;
    end
end

%% plot difference
figure
subplot(1,2,1)
plot(original_input)
hold on
plot(converted_input)
hold off
subplot(1,2,2)
plot(output)

config.err_type = 'NMSE';
error = calculateError(converted_input,original_input,config)