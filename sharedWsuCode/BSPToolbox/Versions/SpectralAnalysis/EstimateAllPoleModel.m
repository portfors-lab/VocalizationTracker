function [a, err] = EstimateAllPoleModel(x, ora)
%EstimateAllPoleModel: Finds the coefficient of an all-pole model
% 
%   [a, err] = EstimateAllPoleModel(x, or)        
%
%   x    Input signal    
%   or   Order of the model. Default=2
%
%   This function finds the coefficients of an 
%   all-pole model for a signal x using the autocorrelation
%
%   a    Model coefficients
%   err  Error associated with the estimate
%
%   Example: Estimate the all-pole model of an 
%   intracranial pressure signal (ICP) sampled at 125 Hz.
%
%      load ICP.mat; 
%      [a, err] = EstimateAllPoleModel(ABP, 2);
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.181, 1996
%
%   Version 1.00 JB
%
%   See also SetCovarMatrix, and SetConMatrix.

% =================================================
% Default parameters
% ================================================

if ( nargin <1 | nargin>2)
    help EstimateAllPoleModel;
    break
end;

or  = 2;
if exist('ora') & ~isempty(ora),
    or = ora;
    end;

%===============================================
% Models Parameters
% ================================================
x   = x(:);
N   = length(x);
if or >= length(x), error('Model order too large'), end

    X   = SetConvMatrix(x, or+1);
    Xq  = X(1:N+or-1,1:or);
    a   = [1; -Xq\X(2:N+or,1)];
    err = abs(X(1:N+or,1)'*X*a);



