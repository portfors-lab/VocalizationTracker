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


%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1 | nargin>4,
    help ECGDetectQSRank;
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
R  = round(fi/ft)                           % Decimation factor
x  = decimate(x,R);                         % Decimate by a factor of R
fs = fs/R;                                  % New fs
LD = length(x);                             % Length of data


%=========================================================================
% Perform R Detection based on Interbeat
%=========================================================================
[R] = ECGDetectRInterbeat(x, fs, fs, 0);


%=========================================================================
% Detrend Data 
%=========================================================================
yntr  = highpass(x, fs, 5, 2, 0);         % Elliptic highpass filter

%=========================================================================
% Detection based on Ranks - QS
%========================================================================
QS = [];                                    % QS wave holder
wl= 10*fs;                                  % Window Size (10 seconds)

for k=1:wl:LD                               
    if ((k+wl-1 < LD-wl))
        yntrs= highpass(yntr(k:k+wl-1), fs, 5, 2, 0);
        thrr  = prctile(yntrs, 10);
        QSs   = DetectMinima(yntrs, thrr);
        QSs   = QSs(:);   
        QS    = [QS; QSs+k-1];
    else
        yntrs = yntr(k:LD);
        thrr  = prctile(yntrs, 10);
        QSs   = DetectMinima(yntrs, thrr);
        QSs   = QSs(:);   
        QS    = [QS; QSs+k-1];
    end;
end;

QS = unique(QS);
QS = sort(QS);

%=========================================================================
% Classification - Q
%=========================================================================
Qs = [];
Q  = [];
for k=1:length(R),
    Qs      = QS(QS < R(k));
    [m, iq] = min(abs(R(k)-Qs));
    Q       = [Q;Qs(iq)];
    Qs      = [];
end;

%=========================================================================
% Classification - S
%=========================================================================
Ss = [];
S  = [];
for k=1:length(R),
    Ss      = QS(QS > R(k));
    [m, iq] = min(abs(R(k)-Ss));
    S       = [S;Ss(iq)];
    Ss      = [];
end;

%=========================================================================
% Correct the effect of the Filtering
%=========================================================================
[min] = DetectMinima(x);

%-------------------------------------------------------------------------
% Correct the effect of the Filtering - Q (adjust time location)
%-------------------------------------------------------------------------
Qi = Q;
Q  = [];
for k=1:length(Qi),
    [m, iq] = min(abs(Qi(k)-min));
    Q       = [Q;min(iq)];
end;

%-------------------------------------------------------------------------
% Correct the effect of the Filtering - S (adjust time location)
%-------------------------------------------------------------------------
Si = S;
S  = [];
for k=1:length(Si),
    [m, iq] = min(abs(Si(k)-min));
    S       = [S;min(iq)];
end;


%=========================================================================
% Interbeat Intervals (Q-Q), (S-S)
%=========================================================================
Q    = unique(Q);                           % Eliminate repetions
Q    = sort(Q);                             % Sort detection
Q2Q  = diff(Q);                           % Calculate Q-Q interval
S    = unique(S);                           % Eliminate repetions
S    = sort(S);                             % Sort detection
S2S  = diff(S);                           % Calculate S-S interval
Q2Qi = (Q(1:length(Q)-1) + Q(2:length(Q)))/2;
S2Si = (S(1:length(S)-1) + S(2:length(S)))/2;


%=========================================================================
% Plotting
%=========================================================================
if pf == 1    
    tx   = 1:length(x);
    
    figure (1)
    figureset(1)
    AxiSset(10);
    h=plot(tx, x, Q, x(Q), 'r.', S, x(S), 'k.');
    set(h, 'MarkerSize', 12);
    axis('tight');
    title('ECG QS-wave Detection: Rank method');
    xlabel('Sample');
    ylabel('Magnitude');
    
    figure (2)
    figureset(2)
    subplot(2,1,1),
    AxiSset(10);
    h=plot(Q2Qi, Q2Q);
    set(h, 'Linewidth', 1.5);
    axis('tight');
    title('ECG QS-wave Detection: Q-Q & S-S Intervals'),
    ylabel('Q-Q Interval');
    
    subplot(2,1,2),
    AxiSset(10);
    h=plot(S2Si, S2S);
    set(h, 'Linewidth', 1.5);
    axis('tight');
    ylabel('S-S Interval');
    xlabel('Sample');
end;


%=====================================================================
% Take care of outputs
%=====================================================================
if nargout==0,
    clear('Q', 'S');
end;   
