% This function is a translation from Spanish to English from BSA algorithm
% (BSA, a Fast and Accurate Spike Train Encoding Scheme) implemented by Imanol
% Bilbao (2018) and translate by Dr. Israel Espinosa
% (https://kedri.aut.ac.nz/staff/staff-profiles/israel-espinosa-ramos)
% Auckland University of Technology, Auckland, New Zealand
% Knowledge Engineering an Discovery Research Institute
% https://kedri.aut.ac.nz/
function [encoded_data, input_mul_data, filter, output_mul_data]= encodeBSAspike(input_signal, filter, output_signal)


    % multiplex signal
    input_mul_data = zeros((size(input_signal,1))*filter.time_period,size(input_signal,2));
    for p = 0:filter.input_length-1
        input_mul_data(mod(1:size(input_mul_data,1),filter.time_period) == p,:) = input_signal;
    end
    
    if nargin > 2
    output_mul_data = zeros((size(output_signal,1))*filter.time_period,size(output_signal,2));
    for p = 0:filter.input_length-1
        output_mul_data(mod(1:size(output_mul_data,1),filter.time_period) == p,:) = output_signal;
    end
    else
        output_mul_data = [];
    end
    
for num_signals = 1:size(input_signal,2)
    
    % Creation of a passband filter
    threshold = filter.threshold;
    filter.f=fir1(filter.order, filter.passband);
    %f=fir1(filter.order, filter.passband)*2;
    %f=fir1(filter.order, filter.passband)*max(data)*2; % Suggestion
    filterSize=length(filter.f);
    
    filter_data=cat(1,(ones(filterSize,1)*input_mul_data(1,num_signals)),input_mul_data(:,num_signals),(ones(filterSize,1)*input_mul_data(end,num_signals))); % add two vectors in the begining and in the end with n-order elements with the last value of the signal
    
    % BSA: The signal must be a vector of n elements
    n_filter = length(filter.f);
    n_input = length(filter_data);
    
    %encoded_data = zeros(,size(input_signal,2))
    for i = 1:n_input
        error1 = 0;
        error2 = 0;
        
        for j = 1:n_filter
            if i+j-1 <= n_input
                error1 = error1 + abs(filter_data(i+j-1) - filter.f(j));
                error2 = error2 + abs(filter_data(i+j-1));
            end
        end
        
        if error1 <= (error2 - threshold)
            %if error1 <= (error2 * threshold) % multiplicative
            encoded_signal(i)=1;
            for j = 1:n_filter
                if i+j-1 <= n_input
                    filter_data(i+j-1) = filter_data(i+j-1)-filter.f(j);
                end
            end
        else
            encoded_signal(i)=0;
        end
    end
    encoded_data(:,num_signals)=encoded_signal';
end

end