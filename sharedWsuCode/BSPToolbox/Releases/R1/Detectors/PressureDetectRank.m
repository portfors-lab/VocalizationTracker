function [B,A,C,D] = PressureDetectRank(x,fsa,dta,pfa);
%PressureDetectRank: Pressure detector based on ranks
%
%   [B,A,C,D] = PressureDetectRank(x,fs,dt,pf);
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=125 Hz      
%   dt   Type: 0=perform detection on raw data, 1=detrend first (default) 
%   pf   Plot flag: 0=none (default), 1=screen,
%
%   B    Percusion peak (P1) index, samples
%   A    Minima preceding the percusion peak (P1) index, samples 
%   C    Dicrotic notch (DN) index, samples
%   D    Dicrotic peak (P2) index, samples
%
%   Beat detection algorithm for pressure signals based on ranks. The 
%   algorithm detects the minima (A) prior to the percursor peak, the 
%   percursor peak (B), the dicrotic notch (C), and the dicrotic peak (D).
%
%   Example: Find the minima (A) prior to the percursor peak, the percursor
%   peak (B), the dicrotic notch (C), and the dicrotic peak (D) in the ICP 
%   signal.
%
%      load ICP; 
%      [B,A,C,D] = PressureDetectRank(icp,fs,1,1);
%
%   Version 1.00 MA
%
%   See also PressureDetect, PressureDetectInterbeat, and ECGDetectQRS.
