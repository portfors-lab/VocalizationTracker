function [Y] = RemoveGaussian(X, wha, wva, dha, dva, pfa)
%RemoveGaussian: Two-Dimensional Gaussian Noise Removing Filter
%
%   [Y] = RemoveGaussian(X,wh,wv,dh,dv,pf)
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
%   Filters a two-dimensional signal using a Mean Filter. The matrix
%   is padded by repeating the values at each edge. Then, a sliding
%   window of specified width and length is placed at the initial value
%   of the input matrix. The mean of the values inside the window is 
%   calculated and stored in the output matrix. The same procedure is 
%   repeated as the window slides through the data. The number of 
%   the window advances at each step in the horizontal and vertical 
%   directions is determined by the values of the input parameters dh
%   and dv, respectively.
%
%   Example: Filter the spectrogram of the ICP signal using a Gaussian
%   Noise Removing Filter with a window of 11-by-21 samples.
%
%       load ICP.mat;
%       icpd = decimate(icp, 15);
%       [S,t,f] = Spectrogram(icpd,125/15,[],[],[],[],1);
%       [Y] = RemoveGaussian(S, 11,21,5,5,1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, 1997.
%
%   Version 1.00 CC
%
%   See also RemoveImpulses. 
