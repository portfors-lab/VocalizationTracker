function [B,t,d,pa] = Beatogram(x,fsa,wla,doa,rga,nsa,pfa);
%Beatogram: Plots pulse morphology versus time as an image
%
%   [B,t,d,pa] = Beatogram(x,fs,do,rg,ns,pf);
%
%   x    Input signal       
%   fs   Sample rate in (Hz). Default = 125 Hz.
%   wl   Length of window to use (s). Default = 8 s. 
%   do   Order of derivative to show (integer 0-2). Default = 0.
%   rg   Range of each beat to show (s). Default = [-1 1].
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        Default = min(2^10,NX).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   B    Matrix containing the image created
%   t    Times (columns) at which B was estimated (s)
%   d    Delays (rows) at which B was estimated (s)
%   pa   Average number of peaks detected in each window
%
%   Generates an image that shows how the pulse morphology changes over
%   time. If the order is specified, then Beatogram plots the 1st or
%   2nd derivative of each pulse.
%
%   Example: Plot the Beatogram of an intra-cranial pressure (ICP)
%   signal.
%
%      load ICP;
%      Beatogram(icp,fs);
%
%   J. McNames, J. Bassale, M. Aboy, C. Crespo, B. Goldstein, "Techniques 
%   for the Visualization of Nonstationary Biomedical Signals," accepted 
%   for presentation at Biosignal 2002. 
%
%   Version 1.00 JM
%
%   See also Spectrogram and Cohereogram.
