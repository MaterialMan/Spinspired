function y = LeakyReLU(x)
%Leaky ReLu node
y = x;
y(x > 0) = x(x > 0);
y(x < 0) = 0.01*x(x < 0);