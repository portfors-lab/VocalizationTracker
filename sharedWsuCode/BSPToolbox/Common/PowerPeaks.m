function [pi,srxsp] = PowerPeaks(x,fsa,fca,pwa,xma,pfa);
%PowerPeaks: Finds peaks in the lowpass filtered signal power 
%
%   [pi,srxsp] = PowerPeaks(x,fs,fc,pw,xmn,pf);
%
%   x      Input signal.
%   fs     Signal sample rate (Hz). Default = 125 Hz.      
%   fc     Cutoff frequency (Hz). Default = fs/4 Hz.
%   pw     Energy power argument. Default = 2.
%   xmn    Minimum value necessary to qualify as an energy peak. 
%          Default = 0.
%   pfa    Plot flag argument: 0=none (default), 1=screen.
% 
%   pi     Vector of signal indices that contain the peaks.
%   srxsp  Smoothed signal power.
%
%   Returns the peaks in the local smoothed signal energy. Signal energy 
%   is defined here as the absolute value of the signal to the power of p. 
%   sorted vector x that are at least as big as xm (optional). Uses a 
%   Lowpass filter type 4 to prevent ringing after impulses.
%
%
%
%   Example: Find the energy peaks in an electrocardiogram signal.
%
%      load ECG;
%      x  = decimate(ecg,10);
%      [pi,y] = PowerPeaks(x,fs/10,2,4);
%      k = 1:length(y);
%      t = (k-1)/fs;
%      plot(t,x,'k',t,y,'b',t(pi),y(pi),'r.','MarkerSize',20);
%      xlim([10 10.5])
%      legend('Signal','Signal Power','Power Peaks');
%
%   Desmond J. Higham and Nicholas J. Higham, "MATLAB Guide," 2000.
%
%   Version 1.00 JM
%
%   See also DetectMaxima, Lowpass, and ECGDetectREnergy.

%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1,
    help PowerPeaks;
    return;
    end;

fs = 125;                                                  % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;
 
fc = fs/4;                                                 % Default cutoff frequency, Hz
if exist('fca') & ~isempty(fca),
    fc = fca;
    end;
    
pw = 4;                                                    % Default signal power
if exist('pwa') & ~isempty(pwa),
    pw = pwa;
end;

xm = 0;                                                    % Default threshold
if exist('xma') & ~isempty(xma),
    xm = xma;
end;    
    
pf = 0;                                                    % Default - no plotting
if nargout==0,                       % Plot if no output arguments
    pf = 1;
end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
end;

%====================================================================
% Core Processing
%====================================================================  
xsp   = Lowpass(abs(x).^pw,fs,fc,4);                       % Smooth the signal.^pw
% xsp   = Lowpass(x.^pw,fs,fc,4);                       % Smooth the signal.^pw
srxsp = xsp.^(1/pw);                                  % Scale down the power of signal
xm = xm./(max(x)/max(srxsp));                         % Scale the threshold to match with plot
pi    = DetectMaxima(srxsp,xm,0);                          % Detect maxima
    
%====================================================================
% Post Processing
%====================================================================  
srxsp = srxsp(:);                                          % Make into a columnvector
srxsp = srxsp.*(max(x)/max(srxsp));                        % Scale srxsp to make it look similar to x

%====================================================================
% Plot Results
%====================================================================  
if pf,
    figure;
    FigureSet;
    k = 1:length(srxsp);
    t = (k-1)/fs;
    xn = srxsp;
    h = plot(t,x,'k',t,xn,'k',t(pi),xn(pi),'k.');
    set(h(1),'Color',0.3*[1 1 1]);                         % Very light gray
    set(h(1),'LineWidth',0.5);
    set(h(1),'Color',0.33*[1 1 1]);                        % Light gray
    set(h(2),'LineWidth',1.5);
    set(h(1),'Color',0.7*[1 1 1]);                         % Gray
    if ~isempty(pi)
        set(h(3),'MarkerSize',20);
    end
    title('Example of Detected Energy Peaks');
    xlabel('Samples');
    ylabel('Signal & Signal Power');
    AxisSet;
    box off;
    legend('Signal','Signal Power','Power Peaks');
    zoom on;
    end;

%====================================================================
% Process Output Arguments
%====================================================================      
if nargout==0,
    clear('pi','srxsp');
    return;
    end;    