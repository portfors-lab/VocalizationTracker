function [A,t,d] = Autocorrelogram(x,SRarg,DRarg,WLarg,NDarg,NSarg,TSarg);
% [A,t,d] = Autocorrelogram(x,fs,DR,WL,ND,NS,TS);
%   Calculates estimates of the spectral content at the specified times
%   using a blackman window of the specified length.
%
%   x : Input signal.
%   fs: Sample rate (default = 1).
%   DR: Array of minimum and maximum delays in seconds (default 0 to (WL-1)/fs).
%   WL: Length of window to use (default 512). Must be even.
%   ND: Number of frequencies to evaluate (default WL)
%   NS: Requested number of times (horizontal pixels) to evaluate 
%       (default 512).
%   TS: Time (in seconds) of the first element of the input signal 
%       (default 0)
