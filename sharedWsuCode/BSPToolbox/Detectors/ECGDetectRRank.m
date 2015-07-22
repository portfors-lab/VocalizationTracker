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
%      load NoisyECG; 
%      [R,nfn,nfp,pr] = ECGDetectRRank(ecg,500,100,1);
%
%   Version 1.00 MA
%
%   See also ECGDetectRInterbeat, ECGDetectR and ECGDetectQRS.


%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1 | nargin>4,
    help ECGDetectRRank;
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


%=========================================================================
% Decimate Input
%=========================================================================
xi = x(:);                                  % Original input, xi
R  = round(fs/ft);                          % Decimation factor
x  = decimate(x,R);                         % Decimate by a factor of R
fs = fs/R;                                  % New fs for detection
LD = length(x);                             % Length of data

%=========================================================================
% Smooth and Detrend Data 
%=========================================================================
yntr  = highpass(x, fs, 5, 2, 0);         % Elliptic highpass filter

%=========================================================================
% Detection based on Ranks
%=========================================================================
R = [];                                     % R wave holder
wl= 5*fs;                                  % Window Size (10 seconds)

for k=1:wl:LD                               
    if ((k+wl-1 < LD-wl))
        yntrs= highpass(yntr(k:k+wl-1), fs, 5, 2, 0);
        thrr = prctile(yntrs, 99.5);
        Rs   = DetectMaxima(yntrs, thrr);
        Rs   = Rs(:);   
        R    = [R; Rs+k-1];
    else
        yntrs= yntr(k:LD);
        thrr = prctile(yntrs, 99.5);
        Rs   = DetectMaxima(yntrs, thrr);
        Rs   = Rs(:);   
        R    = [R; Rs+k-1];
    end;
end;

%=========================================================================
% Interbeat Interval (R-R)
%=========================================================================
R     = unique(R);                           % Eliminate repetions
R     = sort(R);                             % Sort detection
R2R   = diff(R);                             % Calculate R-R interval
R2Ri  = (R(1:length(R)-1) + R(2:length(R)))/2;

%=========================================================================
% Estimated False Positives and False Negatives
%=========================================================================
i     = find((R2R < 0.75*median(R2R)) | (R2R > 1.75*median(R2R)));
nmd   = length(i);
md    = R(i);  
pr    = 1-(nmd/length(R));
ifp   = find(R2R < 0.75*median(R2R));
ifn   = find(R2R > 1.75*median(R2R));
nfp   = length(ifp);
nfn   = length(ifn);


%=========================================================================
% Plotting
%=========================================================================
if pf == 1    
    tx   = 1:length(x);
    
    figure (1)
    figureset(1)
    AxiSset(10);
    h=plot(tx, x, R, x(R), 'r.');
    set(h, 'MarkerSize', 15);
    axis('tight');
    title('ECG R-wave Detection: Rank method');
    xlabel('Sample');
    ylabel('Magnitude');
    box off;
    
    figure (2)
    figureset(2)
    AxiSset(10);
    h=plot(R2Ri, R2R);
    set(h, 'Linewidth', 1.5);
    axis('tight');
    title('ECG R-wave Detection: R-R Intervals');
    xlabel('Sample');
    ylabel('R-R Interval');
    box off;
end;

%=====================================================================
% Take care of outputs
%=====================================================================
if nargout==0,
    clear('R');
end;   
