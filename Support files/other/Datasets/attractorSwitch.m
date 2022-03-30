
function [out] = attractorSwitch(data_length,num_attractors,to_plot,dimension)

rng(1,'twister')
out = [];
t = linspace(-10,10,data_length/num_attractors);
phase = 1;
switch(dimension)
    case 1
        out(:,1) = sin(t*((2*pi)*phase)) + 2*rand*cos(t*((2*pi)*phase))-1;
    case 2
        out(:,1) = sin(t*((2*pi)*phase)) + 2*rand*cos(t*((2*pi)*phase))-1;
        out(:,2) = cos(t*((2*pi)*phase)) + 2*rand*sin(t*((2*pi)*phase))-1;
    case 3
        out(:,1) = sin(t*((2*pi)*phase)) + 2*rand*cos(t*((2*pi)*phase))-1;
        out(:,2) = cos(t*((2*pi)*phase)) + 2*rand*sin(t*((2*pi)*phase))-1;
        out(:,3) = cos(t*((2*pi)*phase)) + 2*rand*sin(t*((2*pi)*phase))-1;
    otherwise
end

for i = 1:num_attractors-1
    phase = 2*rand-1;
    amplitude = 2*rand;
    switch(dimension)
        case 1
            ex = [amplitude*sin(t*(2*pi)*phase) + randi([-8 8])]';
        case 2
            ex = [amplitude*sin(t*((2*pi)*phase)) + randi([-8 8]);...
                amplitude*cos(t*((2*pi)*phase)) + randi([-8 8])]';
        case 3
            ex = [amplitude*sin(t*((2*pi)*phase)) + randi([-8 8]);...
                amplitude*cos(t*((2*pi)*phase)) + randi([-8 8]);...
                amplitude*cos(t*((2*pi)*phase)) + randi([-8 8])]';
        otherwise
    end
    out = [out; ex];
end

if to_plot
    tail = 50;
    for k = tail+1:10:length(out)-1
        switch(dimension)
            case 1
                plot(k-tail:k,out(k-tail:k,1))
                xlim([-10 10])
                ylim([-10 10])
            case 2
                plot(out(k-tail:k,1),out(k-tail:k,2))
                xlim([-10 10])
                ylim([-10 10])
            case 3
                plot3(out(k-tail:k,1),out(k-tail:k,2),out(k-tail:k,3))
                xlim([-10 10])
                ylim([-10 10])
                zlim([-10 10])
            otherwise
                
        end
        drawnow
    end
end