function [R] = ECGDetectRDetector(x,fsa,fsda,pfa);
%ECGDetectRDetector: ECG R-wave detector based on ranks and interbeats
%
%   [R] = ECGDetectRDetector(x,fs,fsd,pf);
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
%
%      load NoisyECG; 
%      [R] = ECGDetectRDetector(ecg,500,500,1);
%
%   Version 1.00 MA
%
%   See also ECGDetectRRank, ECGDetectR and ECGDetectQRS.


%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1 | nargin>4,
    help ECGDetectRInterbeat;
    return;
    end;

fs = 500;                                   % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;
    
fsd = 100;                                  % Detection sampling rate, Hz
ft  = fsd;                                  % Target frequency  
if exist('fsda') & ~isempty(fsda),
    ft = fsda;
    end;
    
pf = 0;                                     % Default - no plotting
if nargout==0,                              % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%===========================================================================
% Process function arguments
%===========================================================================
Lx = length(x);
wl = fs*60*5;
ws = fs*60*1;
f  = [];
for k=1:wl:Lx                              
    if ((k < Lx-wl))
        xs       = x(k:k+wl-1);
        fss = ECGDetectRInterbeat(xs,128,128,0);
        f = [f;fss+k-1];
    else
        xs       = x(k-ws:Lx);
        fss = ECGDetectRInterbeat(xs,128,128,0); 
        f = [f;fss+k-ws-1];
     end;
 end;
f = unique(f);
R = f;
 %================================================
% Plotting
%================================================
if pf == 1,
    t  = (1:length(x))./fs;
    tf = f./fs;
    figure;
    figureset(1, 'wide');
    h1 = axes('Position', [0.10 0.32 0.80 0.60]);     
    h  = plot(t,x,'b',tf,x(f),'r.');
    set(h1, 'Xtick', []);
    set(h, 'Markersize', 15);
    ylabel('Pressure Signal & Detection');
    axis tight;
    axisset(8);
    box off;
    
    h2   = axes('Position', [0.10 0.10 0.80 0.20]);
    f2f  = (diff(f)-1)/fs;
    f2fi = (f(1:length(f)-1)+f(2:length(f)))/2;
    h    = plot(f2fi, f2f, f2fi, f2f, 'b');
    ylabel('IBI,s');
    xlabel('Time, s');
    axis tight;
    axisset(8);
    box off;
end;

%=====================================================================
% Take care of outputs
%=====================================================================
if nargout==0,
    clear('f');
end; 
