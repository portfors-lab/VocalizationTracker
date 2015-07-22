function [C,t,f] = Cohereogram(x1,x2,fsa,sla,wla,sna,ola,fra,nfa,nsa,pfa);
%Cohereogram: Nonstationary estimate the coherency versus time
%
%   [C,t,f] = Cohereogram(x1,x2,fs,sl,wl,sn,ol,fr,nf,ns,pf);
%   
%   x1   First input signal.
%   x2   Second input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   sl   Length of signal segments used to generate estimate (sec). 
%        Default = (signal duration)/100.
%   wl   Length of each of the windows applied to subsets of each 
%        segment (s). Default: sl/20. If vector, wl is used as
%        the window. Otherwise a Blackman window is applied.
%   sn   Signal to noise ratio used to bias the coherence. 
%        Default = inf.
%   ol   Overlap of windows applied to each subset (%). Default = 50%.
%   fr   Minimum and maximum frequencies to display (Hz).
%        Default = [0 fs/2].
%   nf   Number of frequencies to evaluate (vertical resolution). 
%        Default = max(128,round(wl/2)).
%   ns   Requested number of times (horizontal pixels) to evaluate 
%        Default = min(400,length(x)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   C    Matrix containing the image of the cohereogram.
%   t    Times at which the cohereogram was evaluated (s).
%   f    Frequencies at which the cohereogram was evaluated (Hz).
%
%   This function estimates the coherence (not the coherency 
%   spectrum) of two input signals in a moving window and plots the 
%   result as an image.  The coherence is a measure of correlation of 
%   estimated spectral content of two signals as a function of 
%   frequency. The range of values is 0 to 1. 
%
%   Zero padding is used to improve frequency resolution. The 
%   coherence is calculated by estimating the following ratio.
%
%        C = abs(Pxy)/sqrt(Pxx.*Pyy).
%
%   If no output arguments are specified, the image is generated with 
%   a plot of the signals in the bottom axis. The color map has a 
%   fixed range of 0 to 1. The signal is extrapolated at the edges by 
%   repeating the value of the edge.
%
%   Example:  Plot the coherence of 20 minute segments of ABP and ICP
%   data. Limit the frequency axis to only show values between 0 and 
%   4 Hz. Use a segment length of 1 minute and 10 second subsegments.
%   Use a triangular window.
%
%      load ABPICP.mat
%      x1 = decimate(abp,15);
%      x2 = decimate(icp,15);
%      fs = fs/15;
%      wn = triang(round(10*fs));
%      Cohereogram(x1,x2,fs,60,wn);
%
%   Version 1.00 JM
%
%   See also Coherency, Spectrogram, and cohere.
