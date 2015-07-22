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
%      load NoisyECG; 
%      [R,nfn,nfp,pr] = ECGDetectRInterbeat(ecg,500,100,1);
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
    
%=========================================================================
% Decimate Input
%=========================================================================
xi = x(:);                                  % Original input, xi
fi = fs;                                    % Initial frequency
R  = round(fi/ft);                          % Decimation factor
x  = decimate(x,R);                         % Decimate by a factor of R
fs = fs/R;                                  % New fs for detection
LD = length(x);                             % Length of data


%=========================================================================
% Perform R Detection based on ranks
%=========================================================================
[R] = ECGDetectRRank(x, fs, fs, 0);


%=========================================================================
% Interbeat Interval (R-R)
%=========================================================================
R    = unique(R);                           % Eliminate repetions
R    = sort(R);                             % Sort detection
R2R  = diff(R);                             % Calculate R-R interval
R2Ri = (R(1:length(R)-1) + R(2:length(R)))/2;


%=========================================================================
% Correct False Negatives
%=========================================================================
ibi   = median(R2R);
ifn   = find(R2R > 1.5*median(R2R));

mal   = 30*length(ifn);
cnt   = 1;
while ~isempty(ifn)
    ifnp  = ifn(1);
    ni    = R(ifn(1))+1;
    nf    = R(ifn(1)+1)-1;
    fns   = x(ni:nf);
    Rs    = Detectmaxima(fns);
    [m,i] = max(fns(Rs));
    Rs    = ni-1+Rs(i);
    R     = [R; Rs];
    R     = unique(R);
    R     = sort(R);
    R2R   = diff(R);
    ibi   = median(R2R);
    ifn   = find(R2R > 1.5*median(R2R));
    if cnt >= mal,
        break;
    end;
    cnt   = cnt +1;    
end;

%=========================================================================
% Correct False Positives
%=========================================================================
ibi   = median(R2R);
ifp   = find(R2R < 0.75*ibi);
ifp   = ifp+1;

mal   = 30*length(ifp);
cnt   = 1;
while ~isempty(ifp),
    fp    = R(ifp(1));
    R     = R(R~=fp);
    R2R   = diff(R);
    ibi   = median(R2R);
    ifp   = find(R2R < 0.5*ibi);
    ifp   = ifp+1;
    if cnt >= mal,
        break;
    end;
    cnt   = cnt +1;  
end;

%=========================================================================
% Interbeat Interval (R-R)
%=========================================================================
R    = unique(R);                           % Eliminate repetions
R    = sort(R);                             % Sort detection
R2R  = diff(R);                             % Calculate R-R interval
R2Ri = (R(1:length(R)-1) + R(2:length(R)))/2;

%=========================================================================
% Estimated False Positives and False Negatives
%=========================================================================
i    = find((R2R < 0.75*median(R2R)) | (R2R > 1.75*median(R2R)));
nmd  = length(i);
md   = R(i);  
pr   = 1-(nmd/length(R));
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
    axis('tight'); box off; axisset;
    title('ECG R-wave Detection: Interbeat method');
    xlabel('Sample');
    ylabel('Magnitude');
    
    figure (2)
    figureset(2)
    AxiSset(10);
    h=plot(R2Ri, R2R);
    set(h, 'Linewidth', 1.5);
    axis('tight'); box off; axisset;
    title('ECG R-wave Detection: R-R Intervals');
    xlabel('Sample');
    ylabel('R-R Interval');  
end;



%=====================================================================
% Take care of outputs
%=====================================================================
if nargout==0,
    clear('R');
end;   
















