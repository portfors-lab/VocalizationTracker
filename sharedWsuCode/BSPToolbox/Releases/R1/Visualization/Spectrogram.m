function [S,t,f] = Spectrogram(x,fsa,wla,fra,nfa,nsa,pfa);
%Spectrogram: Generates the spectrogram of the specified signal
%
%   [S,t,f] = Spectrogram(x,fs,wl,fr,nf,ns,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1024 samples.
%        If a vector, specifies entire window.
%   fr   Minimum and maximum frequencies to display (Hz).
%        Default = [0 fs/2].
%   nf   Number of frequencies to evaluate (vertical resolution). 
%        Default = max(128,round(wl/2)).
%   ns   Requested number of times (horizontal pixels) to evaluate 
%        Default = min(400,length(x)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   S    Matrix containing the image of the Spectrogram.
%   t    Times at which the spectrogram was evaluated (s).
%   f    Frequencies at which the spectrogram was evaluated (Hz).
%
%   Calculates estimates of the spectral content at the specified 
%   times using a modified periodogram. The mean of the signal is 
%   removed as a preprocessing step. The square root of the power 
%   spectral density is calculated and displayed. To limit 
%   computation, decimate the signal if necessary to make the upper 
%   frequency range approximately equal to half the sampling 
%   frequency. 
%   
%   If only the window length is specified, the blackman window is 
%   used. The specified window length should be odd. If called with
%   a window with an even number of elements, a zero is appended to
%   make the window odd.
%
%   Example: Generate the spectrogram of an intracranial pressure
%   signal using a Blackman-Harris window that is 45 s in duration.
%
%      load ICP.mat; 
%      icpd = decimate(icp,15);
%      wl   = round(45*fs/15);
%      Spectrogram(icpd,fs/15,blackmanharris(wl));
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, 1996.
%
%   Version 1.00 JM
%
%   See also specgram, window, decimate, Beatogram, and Cohereogram.
