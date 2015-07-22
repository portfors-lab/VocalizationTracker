function [B,A,C,D] = PressureDetectInterbeat(x,fsa,wla,pfa);
%PressureDetect: Pressure detector algorithm (ICP & ABP)
%
%   [B, A, C, D] = PressureDetect(x,fs,wl,pf);
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=125 Hz      
%   wl   Window length (s), default=5 s 
%   pf   Plot flag: 0=none (default), 1=screen
%
%   B    Percusion peak (P1) index, samples
%   A    Minima preceding the percusion peak (P1) index, samples
%   C    Dichrotic notch (DN) index, samples
%   D    Dichrotic peak (P2) index, samples
%
%   Beat detection algorithm for pressure signals based on ranks and 
%   intebeat distance. The algorithm detects the minima prior to the 
%   percursor peak (A), the percursor peak (B), the dicrotic notch (C),
%   and the dicrotic peak (D). The algorithm (1) detects all maxima and
%   minima in the raw data, (2) detrends and smooths the data using a 
%   bandpass filter, (3) estimates the heart rate of the input signal as
%   a function of time (window), (4) detects and classifies all maxima 
%   and minima in each window of filtered signal based on ranks, and 
%   (5) corrects the missed beats and false positives based on the 
%   interbeat distance. 
%
%   Example: Find the minima (A) prior to the percursor peak, the percursor
%   peak (B), the dicrotic notch (C), and the dicrotic peak (D) in the ICP 
%   signal
%
%      load ICPCOMPOSITE; 
%      PressureDetect(icpcomposite, fs);
%
%   Version 1.0 MA
%
%   See also PressureDetectInterbeat, PressureDetectRank, and ECGDetectQRS.
