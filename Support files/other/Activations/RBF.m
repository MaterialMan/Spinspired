function Y = RBF(X,Beta,w)

if nargin < 2
    Beta = 1;
    w = 0;
end

Y = exp(-Beta.*norm(X - w).^2);

