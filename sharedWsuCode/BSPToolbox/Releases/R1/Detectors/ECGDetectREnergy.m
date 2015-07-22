function [t,r] = ECGDetectREnergy(x,fsa,hra,pfa);
%DetectQRSEnergyBased: Detects QRS components of ECG signals
%
%   [t,a] = ECGDetectqREnergy(x,fs,hr,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=500 Hz.
%   hr   Expected range of the heart rate (Hz). Default=[0.8 3.4] Hz.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   t    Matrix of times (sec) of R components
%   r    Matrix of amplitudes  of R components 
%
%   Identifies the R time and amplitude for each individual beat. 
%   The amplitudes are measured relative to the detrended signal.
%   This is a complicated algorithm that attempts to locate the QRS
%   complexes based on amplitude and then corrects the times of the
%   estimated complexes using inter-beat variability.
%
%   Note that the times and amplitudes are returned because the
%   signal is internally resampled at approximately 1000 Hz. 
%
%   Example: Estimate ECG component times and amplitudes of an 
%   ECG signal.
%
%      load NOISYECG.mat; 
%      [t,a] = ECGDetectREnergy(noisyecg,fs,[],1);
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.408-412, 1996.
%
%   Version 1.00 JM
%
%   See also Detectors and HarmonicPSD.
