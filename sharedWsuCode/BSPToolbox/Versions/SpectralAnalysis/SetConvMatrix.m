function [M] = SetConvMatrix(x,ora)
%SetConvMatrix: This function sets up a convolution matrix
% 
%   [M] = SetConvMatrix(x,or);       
%
%   x    Input signal    
%   or   Order of the model. Default=2
%
%   This function sets up a convolution matrix
%  
%   Example: Sets a 2 order convolution matrix
%
%      load ICP.mat; 
%      [M] = SetConvMatrix(icp,2);
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.573, 1996.
%
%   Version 1.00 JB
%
%   See also SetCovarMatrix.


% ==================================================================
% Default parameters
% ==================================================================
if ( nargin < 1 | nargin > 2)
    help SetConvMatrix;
    return;
end;

or  = 2;
if exist('ora') & ~isempty(ora),
    or = ora;
    end;

% ==================================================================
% Estimate Model parameters
% ==================================================================

N = length(x) + 2*or-2;
x = x(:);
xpad = [zeros(or-1,1); x; zeros(or-1,1)];
for i = 1:or
    M(:,i) = xpad(or-i+1:N-i+1);
end;

