function [R,nfn,nfp,pr] = ECGDetectRInterbeat(x,fsa,fsda,pfa);
%ECGDetectRInterbeat: ECG R-wave detector based on ranks and interbeats
%
%   [R,nfn,nfp,pr] = ECGDetectRInterbeat(x,fs,fsd,pf);
%
%   x     Input signal       
%   fs    Signal sampling frequency (Hz). Default=500 Hz      
%   fsd   Signal sampling frequency during detection. Default=100 Hz 
%   pf    Plot flag: 0=none (default), 1=screen
%
%   R     R wave index, samples
%   nfn   Estimated number of false negatives (missed detections)
%   nfp   Estimated number of false positives (over detections)
%   pr    Estimated performance 
%
%   Beat detection algorithm for ECG signals based on interbeats. The 
%   algorithm detects the time location of the R wave based on ranks, 
%   and the interbeat intervals. In order to perform the detection 
%   faster, the algorithm decimates the input signal to fsd, which 
%   is 100 Hz by default. If time resolution of the detected R-wave 
%   is more importand than speed, then fsd can be set to fs, and the 
%   algorithm would not decimate the signal during detection. The 
%   algorithm also outputs an estimate of the number of false negatives, 
%   false positives, and the performance based on the interbeat intervals. 
%
%   Example: Detect the R component in an ECG signal sampled at 500 Hz.
%   Perform the detection at fsd = 100 Hz (the signal is decimated to 
%   have a sampling rate of 100 Hz during detection for fater detection).
%
%      load NOISYECG; 
%      [R,nfn,nfp,pr] = ECGDetectRInterbeat(noisyecg,500,100,1);
%
%   Version 1.00 MA
%
%   See also ECGDetectRRank, ECGDetectR and ECGDetectQRS.
