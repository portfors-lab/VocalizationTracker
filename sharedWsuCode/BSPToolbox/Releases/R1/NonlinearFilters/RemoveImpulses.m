function [Y] = RemoveImpulses(X, wha, wva, dha, dva, pfa)
%RemoveImpulses: Two-Dimensional Impulse Removing Filter
%
%   [Y] = RemoveImpulses(X,wh,wv,dh,dv,pf)
%
%   X       Input signal(must be a matrix)
%   wh      Horizontal size of the sliding window in samples (odd 
%           integer). Default=21.
%   wv      Vertical size of the sliding window in samples (odd 
%           integer). Default=21.
%   dh      Size of the horizontal step of the sliding window (samples).
%           Default=1.
%   dv      Size of the vertical step of the sliding window (samples). 
%           Default=1.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   Y       Filtered signal
%
%   Filters a two-dimensional signal using a Median Filter. The matrix
%   is padded by repeating the values at each edge. Then, a sliding
%   window of specified width and length is placed at the initial value
%   of the input matrix. The median of the values inside the window is 
%   calculated and stored in the output matrix. The same procedure is 
%   repeated as the window slides through the data. The number of 
%   the window advances at each step in the horizontal and vertical 
%   directions is determined by the values of the input parameters dh
%   and dv, respectively.
%
%   Example: Filter the spectrogram of the ICP signal using an Impulse
%   Removing Filter with a window of 11-by-81 samples.
%
%       load ICP.mat;
%       icpd = decimate(icp, 15);
%       [S,t,f] = Spectrogram(icpd,125/15,[],[],[],[],1);
%       [Y] = RemoveImpulses(S, 11,81,10,10,1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, 1997.
%
%   Version 1.00 CC
%
%   See also MedianFilter and RemoveGaussian.
