function [input,output] = generateHenonMap(dataLength, stdev)

i=1;

while(1)
    rng(i,'twister');
    
    noise = stdev*randn(dataLength,1);
    %noise = stdev*rand(dataLength,1);% taken from "A Comparative Study of Reservoir Computing...
    ...for Temporal Signal Processing (Goudarzi,2015)"
 
%x = zeros(dataLength,1);
y = zeros(dataLength+3,1);


for i = 3:dataLength+3
    %y(i) = 1-(1.4*(y(i-1).^2)) + (0.3*y(i-2));
    
    y(i+1) = 1-((1.4)*(y(i).^2)) + (0.3*y(i-1));% + (stdev*randn);
    
    %x(i+1) = (0.3*x(i-1) + (stdev*randn)) - 1.4*(x(i).^2);
    
    % x(i+1)= 1-1.4*x(i).^2 + y(i);
    % y(i+1)=0.3*x(i);%+ noise(i+1);
    %y(i+1)= (y(i+1)-0.5)*2;
    
    %x(i+1)= y(i) - 1.4*(x(i).^2);
    %y(i+1)= 0.3*x(i) + stdev*randn;
    
end

%y = ((y+noise)-0.5)*2;
%y = y+noise;
%input = [0; y(1:dataLength-1)];
%output = y;

y = y(4:end);

ahead = 1;
input = y(1:end-ahead);
output = y(ahead+1:end);
        
if isinf(output(end-1))
    i = i+1;
else
    break
end
end

% figure
% scatter(input,output)