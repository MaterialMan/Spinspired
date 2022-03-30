
function y = mackeyGlassNode(x)

%y = 0.4 * (x./(1+x.^1));
b = 0.4; p = 6.88; C = 1.33;

y = C * (x./(1+b^p *(x.^p)));

% b = 0.4;
% p = 6.88;
% 
% T = 1;
% tau = 80;
% C = 1.33;
% alpha = 1; %gain
% beta = C*alpha;
% 
% y = zeros(length(x),1);
% 
% for t = tau+1:length(x)
%     
%     f_t = C * (alpha*y(t-tau) + beta*x(t));
%     f_b = (1 + b^p) * (alpha*y(t-tau) + beta*x(t))^p;
%     f = f_t/f_b;
%     y(t) = (1/T)*(-y(t) + f);
%     
% end

