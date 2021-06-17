function y = step(x)

for i = 1:length(x)
if mod(floor(x(i)),2) == 0 
    y(i) = 1;
else
    y(i) = -1;
end
end