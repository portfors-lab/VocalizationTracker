function [hp,fh] = HarmonicPSD(p,fa,nha,tha,pfa)
%HarmonicPSD: Estimate the total (harmonic) PSD versus frequency.
%
%   [hp,fh] = HarmonicPSD(p,f,th,pf)
%
%   p    Estimated power spectral density (PSD)  
%   f    Vector of frequencies at which p was estimated at. 
%        Default = 0 to almost 1.0 Hz (equivalent fs = 2 Hz).
%   nh   Number of harmonics to add. Default = 5.
%   th   Threshold for fraction of harmonic power added. Default = 2.     
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   hp   Estimated harmonic PSD (column vector).
%   fh   Frequencies at which estimate was made.
%
%   Estimates the total power at each frequency including the power
%   of the harmonics. This method was described in detail in the
%   reference below.
%
%   Example: Estimate the harmonic PSD of an ECG signal signal 
%   sampled at 500 Hz using 11 harmonics and plot the results.
%
%      load ECG;
%      [p,f] = BlackmanTukey(ecg-mean(ecg),fs,20);
%      HarmonicPSD(p,f,11);
%
%   J. McNames, C. Crespo, M. Aboy, J. Bassale, L. Jenkins, 
%   B. Goldstein, "Harmonic Spectrogram for the Analysis of 
%   Semi-periodic Physiologic Signals," submitted to 
%   EMBS and BME 2002. 
%
%   Version 1.00 JM
%
%   See also SPECTRUM, WINDOW, and SpectralAnalysis.
