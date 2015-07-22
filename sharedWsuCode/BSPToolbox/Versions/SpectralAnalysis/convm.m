function X = convm(x,pa)
% This function sets up a convolution matrix
%
%  [X]       : convm(x,p);
%  x         : input signal       
%  p         : order of the model, default 2
%
% This function sets up a convolution matrix
%  
% Examples
%   n = 0:511;
%   T = 1/32;
%   fs = 32;
%   x                   = 0.01*randn(1, 512) + sin(2*pi*2*n*T);
%   [X]                 = convm(x,2);

%  Version 1.02

if ( nargin < 1 | nargin > 2)
    help convm;
    return;
end;

p  = 2;
if exist('pa') & ~isempty(pa),
    p = pa;
    end;


N = length(x) + 2*p-2;
x = x(:);
xpad = [zeros(p-1,1); x; zeros(p-1,1)];
for i = 1:p
    X(:,i) = xpad(p-i+1:N-i+1);
end;

