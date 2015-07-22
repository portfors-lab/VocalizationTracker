function [R,nfn,nfp,pr] = ECGDetectRRank(x,fsa,fsda,pfa);
%ECGDetectRRank: ECG R-wave detector based on ranks
%
%   [R,nfn,nfp,pr] = ECGDetectRRank(x,fs,fsd,pf);
%
%   x     Input signal       
%   fs    Sampling frequency (Hz). Default=500 Hz      
%   fsd   Sampling frequency during detection. Default = 100 Hz 
%   pf    Plot flag: 0=none (default), 1=screen
%
%   R     R wave index, samples
%   nfn   Estimated number of false negatives (missed detections)
%   nfp   Estimated number of false positives (over detections)
%   pr    Estimated performance 
%
%   Beat detection algorithm for ECG signals based on ranks. The 
%   algorithm detects the time location of the R wave. In order to 
%   perform the detection faster, the algorithm decimates the input 
%   signal to fsd, which is 100 Hz by default. The algorithm also 
%   outputs an estimate of the number of false negatives, false 
%   positives, and the performance based on the interbeat intervals. 
%
%   Example: Detect the R component on an ECG signal sampled at
%   500 Hz. Perform the detection at fsd = 100 Hz (the signal is 
%   decimated to have a sampling rate of 100 Hz during detection 
%   for fater detection).
%
%      load NOISYECG; 
%      [R,nfn,nfp,pr] = ECGDetectRRank(noisyecg,500,100,1);
%
%   Version 1.00 MA
%
%   See also ECGDetectRInterbeat, ECGDetectR and ECGDetectQRS.
