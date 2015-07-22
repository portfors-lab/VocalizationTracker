function [a, err] = acm(x,pa)
% This function finds the coefficient of an all-pole model for a signal x(n) using the autocorrelation
%  
%  [a, err]  : acm(x,p)
%
%  x         : input signal       
%  p         : order of the model, default 2
%
%  This function finds the coefficients of an all-pole model for a signal x(n) using the autocorrelation
%
%  Example: 
%   n = 0:511;
%   T = 1/32;
%   fs = 32;
%   x                   = 0.01*randn(1, 512) + sin(2*pi*2*n*T);
%   [a, e]              = acm(x,2);

%  Version 1.02


if ( nargin <1 | nargin>2)
    help acm;
    return;
end;

p  = 2;
if exist('pa') & ~isempty(pa),
    p = pa;
    end;


x   = x(:);
N   = length(x);
if p >= length(x), error('Model order too large'), end
    X   = convm(x, p+1);
    Xq  = X(1:N+p-1,1:p);
    a   = [1; -Xq\X(2:N+p,1)];
    err = abs(X(1:N+p,1)'*X*a);
end


