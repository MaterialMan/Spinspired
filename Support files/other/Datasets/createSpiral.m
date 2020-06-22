function [xy, class] = createSpiral(numSamples,noise)

[x1, y1] = genSpiral(0,numSamples,noise); % Positive examples.
class1 = ones(numSamples/2,1);
[x2, y2] = genSpiral(pi,numSamples,noise); % Negative examples.
class2 = ones(numSamples/2,1).*-1;

xy = [x1 x2; y1 y2]'; 
class = [class1; class2];

    function [x,y] = genSpiral(deltaT, numSamples, noise)
        n = numSamples / 2;
        for i = 1:n
            r = i / n * 5;
            t = 1.75 * i / n * 2 * pi + deltaT;
            x(i) = r * sin(t) + (2*rand-1)*noise;
            y(i) = r * cos(t) + (2*rand-1)*noise;
        end
    end  
end

% 
%     for i = 1:400
%         [x1(i), y1(i)] = getSpiral(i, spiral_num);
%         class(i,1) = spiral_num;
%     end
%     xy = [x1; y1]'; 
%     
%     function [x, y] = getSpiral(i, spiral_num)
%         phi = i/16 * pi;
%         r = 6.5 * ((104 - i)/104);
%         x = (r * cos(phi) * spiral_num)/13 + 0.5;
%         y = (r * sin(phi) * spiral_num)/13 + 0.5;
%     end