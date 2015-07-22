function [Q,S,R] = ECGDetectQSRank(x,fsa,fsda,pfa);
%ECGDetectQSRank: ECG QS-wave detector based on ranks
%
%   [R,s] = ECGDetectQSRank(x,fs,fsd,pf);
%
%   x     Input signal       
%   fs    Signal sample rate (Hz). Default=500 Hz      
%   fsd   Sampling frequency during detection. Default = 100 Hz 
%   pf    Plot flag: 0=none (default), 1=screen,
%
%   Q     Q wave index, samples
%   S     S wave index, samples
%
%   Beat detection algorithm for ECG signals based on ranks. The 
%   algorithm detects the time location of the Q and the S wave. 
%   In order to perform the detection faster, the algorithm decimates
%   the input signal to fsd, which is 100 Hz by default. If time 
%   resolution of the detected Q and S-waves is more importand than 
%   speed, then fsd can be set to fs, and the algorithm would not 
%   decimate the signal during detection. The algorithm also outputs
%   an estimate of the number of false negatives, false positives, 
%   and the performance based on the interbeat intervals. 
%
%   Example: Detect the Q and S waves in an ECG signal sampled at 500 Hz.
%   Perform the detection at fsd = 100 Hz (the signal is decimated to 
%   have a sampling rate of 100 Hz during detection for fater detection).
%
%      load ECG; 
%      [Q,S] = ECGDetectQSRank(ecg,500,100,1);
%
%   Version 1.00 MA
%
%   See also ECGDetectRInterbeat, ECGDetectR and ECGDetectQRS.
