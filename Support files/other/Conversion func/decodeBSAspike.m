function decoded_data = decodeBSAspike(filter,encoded_data,wash_out)

for num_signals = 1:size(encoded_data,2)
    
    data = encoded_data(:,num_signals);
    filterSize=length(filter.f);
    data=conv(data,filter.f);
    data=data(filterSize+1:end-((2*filterSize)-1),:);
    
    if filter.time_period > 1
        data = data(mod(1:size(data,1),filter.time_period) == 0,:);
    end
    
    decoded_data(:,num_signals) = data;
end

decoded_data = decoded_data(wash_out+1:end,:);

%mse=(mean((filter_data(:,:) - decodedData(:,:)).^2)); % mean square error
