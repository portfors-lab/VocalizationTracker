function [C,t,d] = CrossCorrelogram(x1,x2,fsa,dra,wla,nfa,nsa,tsa,pfa);
% CrossCorrelogram: Estimate and plot the cross-correlation of 
% two non-stationary signals.
%
%   [C,t,d] = CrossCorrelogram(x1,x2,fs,dr,wl,nf,ns,ts,pf);
%
%   x1   Input signal.
%   x2   Input signal.
%   fs   Sample rate.  Default=1.
%   dr   Maximum delays in seconds on verticle axis.  Default=wl/2.
%   wl   Length of window to use.  Default=512. 
%   nf   Number of frequencies to evaluate.  Default=wl.
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        default = 1024.
%   ts   Time (in seconds) of the first element of the input signal, 
%        Default=0.
%   pf   Plot flag:  0=none (default), 1=screen.
% 
%   C    CrossCorrelogram output matrix
%   t    Times at which the CrossCorrelogram was evaluated (s).
%   d    Delays at which the CrossCorrelogram was evaluated (s).
%
%   Plots the estimated cross-correlation of two non-stationary 
%   signals, x1 and x2.  The cross-correlation is calculated 
%   across overlapping windows of the signals and plotted using 
%   a surface plot.  The top and bottom 2.5 percentile of each 
%   cross-correlation is removed before plotting.  The x-axis 
%   represents the time in seconds or minutes, depending on the 
%   length of the signal.  The y-axis represents the lag in seconds.
%   The colorbar represents the magnitude of the cross-correlation, 
%   and the two original signals are plotted beneath the correlogram.  
%
%   Example:  Plot the Cross-Correlogram of an ABP and ICP signal 
%   with a sampling rate of 125 Hz and delay window of -1:1.
%
%      load ABPICP.mat
%      [C,t,d] = CrossCorrelogram(abp,icp,125,1,[],[],[],[],1);  
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 1, Time-domain Methods", Med. & Biol. Eng.
%   & Comput., 1990, 28, pp. 509-524.
%
%   Version 1.00 LJ
%
%   See also AutoCorrelogram, Autocorrelation, and CrossCorrelation.
