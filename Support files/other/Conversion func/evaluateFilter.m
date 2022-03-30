function [nmse,encoded_data,decoded_data,f,filter_data] = evaluateFilter(filter,data)

wash_out = 0;

% Encoding the data
[encoded_data,filter_data,filter]= encodeBSAspike(data, filter);

% Signla reconstruction by convolution of the data and the filter. The
% decodedData length dl = el+fl-1 where el is the encodedData length and fl
% the filter length.
decoded_data = decodeBSAspike(filter,encoded_data,wash_out);

if filter.time_period > 1
   filter_data = filter_data(mod(1:size(filter_data,1),filter.time_period) == 0,:);
end
    
nmse= mean((filter_data - decoded_data).^2)/var(filter_data); % mean square error

f = filter.f;
end