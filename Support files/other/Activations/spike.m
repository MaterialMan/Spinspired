function y = spike(x)

for i = 1:length(x)
if mod(floor(x(i)),2) == 0 
    y(i) = 1 - 2*(x(i)-floor(x(i)));
else
    y(i) = -1 + 2*(x(i)-floor(x(i)));
end
end