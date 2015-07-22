function [R] = SetCovarMatrix(x,ora)
%SetCovarMatrix: This function sets up a covariance matrix
% 
%   [R] = SetCovarMatrix(x,or);       
%
%   x    Input signal    
%   or   Order of the model. Default=2
%
%   This function sets up a convariance matrix
%  
%   Example: Sets a 2 order convariance matrix
%
%      load ICP.mat; 
%      [M] = SetCovarMatrix(icp,2);
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.573, 1996.
%
%   Version 1.00 JB
%
%   See also SetConvMatrix.

% ==================================================================
% Default parameters
% ==================================================================

if ( nargin < 1 | nargin > 2)
    help SetCovarMatrix;
    return;
end;

or  = 2;
if exist('ora') & ~isempty(ora),
    or = ora;
    end;

% =================================================================
% Estimate Model parameters
% =================================================================

x = x(:);
m = length(x);

x = x - ones(m,1) * sum((x)/m);
R = SetConvMatrix(x, or+1)'*SetConvMatrix(x, or+1)/(m-1);
