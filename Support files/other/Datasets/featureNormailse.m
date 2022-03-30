function [Y] = featureNormailse(X,config)

%Number of observations
N=length(X(:,1));

%Number of variables
M=length(X(1,:));

% output array of normalised values
Y=zeros(N,M);

switch config.preprocess
    
    case'std'
        %Subtract mean of each Column from data
        Y=X-repmat(mean(X),N,1);
        
        %normalize each observation by the standard deviation of that variable
        Y=Y./repmat(std(X,0,1),N,1);
        
    case 'scaling' % scal each column/feature to be within [0 1]
        
        % determine the maximum value of each colunm of an array
        Max=max(X);
        % determine the minimum value of each colunm of an array
        Min=min(X);
        %array that contains the different between the maximum and minimum value for each column
        Difference=Max-Min;    
        %subtract the minimum value for each column
        Y=X-repmat(Min,N,1);
        %Column by the difference between the maximum and minimum value 
        Y=Y./repmat(Difference,N,1);        
        
    case 'rescale' % all scaled together within [0 1]
        Y =(X-min(min(X)))./(max(max(X))-min(min(X)));
        
    case 'rescale_diff'
        
        % normalise values between -1 and 1
        if abs(min(min(X))) > max(max(X))
            max_range_value = abs(min(min(X)));
            min_range_value = min(min(X));
        else
            max_range_value = max(max(X));
            min_range_value = -max(max(X));
        end
        
        Y = 2 .* X./(max_range_value-min_range_value);       
        
    case 'clip'
        X(X>1) = 1;
        X(X<0) = 0;
        Y = X;
        
    otherwise
        Y = X;
end

% shift data range
if ~isempty(config.preprocess_shift)
    range = config.preprocess_shift;
    Y = (range(1) + (range(2)-range(1))*Y);
end

end