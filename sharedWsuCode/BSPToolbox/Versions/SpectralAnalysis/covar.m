function R = covar(x,pa)
% This function sets up a covariance matrix
%
%  [R]        : convm(x,p);
%   x         : input signal       
%   p         : order of the model, default 2
%
% This function sets up a covariance matrix
%  
% Examples
%   n = 0:511;
%   T = 1/32;
%   fs = 32;
%   x                   = 0.01*randn(1, 512) + sin(2*pi*2*n*T);
%   [X]                 = covar(x,2);

%  Version 1.02

if ( nargin < 1 | nargin > 2)
    help covar;
    return;
end;

p  = 2;
if exist('pa') & ~isempty(pa),
    p = pa;
    end;


x = x(:);
m = length(x);

x = x - ones(m,1) * sum((x)/m);
R = convm(x, p+1)'*convm(x, p+1)/(m-1);
