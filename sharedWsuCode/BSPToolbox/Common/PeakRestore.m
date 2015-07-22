function [y] = PeakRestore(x,fsa,fca,tha,ppa,gfa,pfa)
%PeakRestore: Nonlinear lowpass filter without peak attenuation
%
%   [y] = PeakRestore(x,fs,fc,th,pp,gc,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default=2 Hz.
%   fc   Cutoff frequency (Hz). Default=fs/4.
%   th   Vector containing the two peak thresholds or indices. 
%   pp   Peak polarity. 1=Positive peaks (default), 2=Negative peaks,
%        3 = Both.
%   gf   Gain filter coefficient. Default=0.9.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    Filter output.
%  
%   Linear-phase lowpass filters are often used as a preprocessing 
%   step in signal processing applications to eliminate high 
%   frequency components of noise that do not overlap with the signal 
%   spectrum. Signals with sharp peaks confound this type of noise 
%   reduction because sharp peaks have a broad spectrum and are 
%   significantly attenuated by lowpass filters. This implements 
%   a simple method that restores the signal peaks to their full 
%   amplitude by adding a filtered estimate of the peak residuals to 
%   the original lowpass filter output.
%
%   The fourth input argument is used differently depending on how
%   long the th vector is.
%
%   Scalar:           Any peak that is above th is labeled a peak.
%   2-element vector: Any peak that is between the two elements is
%                     labeled a peak.
%   n-element vector: vector containing the indices of the peaks.
%                     This enables the function to be used with any
%                     peak detector.
%
%   Example: Apply a linear and nonlinear filter to an 
%   electrocardiogram signal.
%
%      load('ECG.mat');   
%      x   = ecg(1:10000);
%      x   = Highpass(x,fs,5); % Drift filter
%      xlp = Lowpass(x,fs,20,4,2);
%      xpr = PeakRestore(x,fs,20);
%      k   = 1:length(x);
%      t   = (k-1)/fs;
%      plot(t,x,'g',t,xlp,'r',t,xpr,'b');
%      xlim([9.4 10.1]);
%      legend('No drift','Linear Lowpassed','Peak Restored',4);
%
%   J. McNames and B. Goldstien, "A Nonlinear Lowpass Filter that 
%   Eliminates Peak Attenuation," ICASSP 2002, in press.
%
%   Version 1.00 JM
% 
%   See also Lowpass.

%====================================================================
% Process function arguments
%====================================================================
if nargin<1,
    help PeakRestore;
    return;
    end;
    
fs = 2;                                                    % Default sampling rate, Hz
if exist('fsa','var') && ~isempty(fsa),
    fs = fsa;
    end;
 
fc = fs/4;                                                 % Default cutoff frequency, Hz
if exist('fca','var') && ~isempty(fca),
    fc = fca;
    end;

tc = 0;                                                    % Threshold code: thresholds are set
th = [prctile(x,95) max(x)];                               % Default (dumb) thresholds
if exist('tha','var') && ~isempty(tha),
    if length(tha)==1
        th = [tha max(x)+std(x)];                          % Automatically pick upper threshold
        tc = 0;                                            % Thresholds code: thresholds are set
    elseif length(tha)==2,
        th = tha;
        tc = 0;                                            % Threshold code: thresholds are set
    else
        pd = th;                                           % User specified peaks 
        tc = 1;                                            % Threshold code: Vector th contains the peaks
        end;
    end;
    
pp = 1;                                                   % Default polarity type
if exist('ppa','var') && ~isempty(ppa),
    pp = ppa;
    end;

gf = 0.9;
if exist('gfa','var') && ~isempty(gfa),
    gf = gfa;
    end;
    
pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa','var') && ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Preprocessing & Memory allocation
%====================================================================  
x  = x(:);                                                 % Convert x to a column vector
nx = length(x);
ft = 4;                                                    % Lowpass Filter type
cf = 2;                                                    % Use noncausal lowpass filter

%====================================================================  
% Main Loop
%====================================================================  
xf = Lowpass(x,fs,fc,ft,cf);                               % Step 1: Filter the signal
e  = x - xf;                                               % Step 2: Calculate residual
ef = Lowpass(e,fs,fc,ft,cf);                               % Step 3: Calculate filtered error

em = zeros(size(ef));                                      % Error masked and amplified

%====================================================================  
% Steps 4 & 5: Calculate Estimated Upper Peak Residuals
%====================================================================  
if pp==1 || pp==3,                                         % Positive Peaks
    if tc==0,                                              % Vector th contains lower and upper threshold
        pi = DetectMaxima(ef);
        pd = pi(ef(pi)>th(1) & ef(pi)<th(2));              % Detected peaks
    elseif tc==1,                                          % User has already specified peaks
        pd = th;                                           % Detected peaks
        end;
    gn = 1;                                                % Initial gain filter value
    for c1 = 1:length(pd),
        id = pd(c1);
        r = e(id)/ef(id);                                  % Estimated Gain Ratio
        if r>0.5,                                          % Only include if gain is greater than 1
            gn = gf*gn + (1-gf)*r;
            end;
        im1 = max(id-1,1);                                 % Find first zero to left
        while im1-1>1 && ef(im1-1)>0,
            im1 = im1 - 1;
            end;
        im2 = min(id+1,nx);                                % Find first zero to right
        while im2+1<nx && ef(im2+1)>0,
            im2 = im2 + 1;
            end;
        em(im1:im2) = em(im1:im2) + gn*e(im1:im2);         % Add filtered & scaled peak
        end;
    end;

%====================================================================  
% Steps 4 & 5: Calculate Estimated Lower Peak Residuals
%====================================================================      
if pp==2 || pp==3,                                         % Negative peaks
    if tc==0,                                              % Vector th contains lower and upper threshold
        pi = DetectMaxima(ef);
        pd = pi(ef(pi)>th(1) & ef(pi)<th(2));              % Detected peaks
    elseif tc==1,                                          % User has already specified peaks
        pd = th;                                           % Detected peaks
        end;
    gn = 1;
    for c1 = 1:length(pd),
        id = pd(c1);
        r = e(id)/ef(id);                                  % Ratio
        if r>0.5,                                          % Only include if gain is greater than 1
            gn = gf*gn + (1-gf)*r;
            end;
        im1 = max(id-1,1);                                 % Find first zero to left
        while im1-1>1 && ef(im1-1)<0,
            im1 = im1 - 1;
            end;
        im2 = min(id+1,nx);                                % Find first zero to right
        while im2+1<nx && ef(im2+1)<0,
            im2 = im2 + 1;
            end;
        em(im1:im2) = em(im1:im2) + gn*e(im1:im2);         % Add filtered & scaled peak
        end;
    end;

%====================================================================  
% Filter the masked residuals & Add
%==================================================================== 
en = Lowpass(em,fs,fc,ft,cf);                              % Step 7: Error nonlinear-filtered   
y  = xf + en;                                              % Step 8: Add nonlinear residual

%==============================================================================
% Plot Example
%==============================================================================
if pf==1,
    figure;
    if exist('FigureSet','file'),
        FigureSet;
        end;
    k = 1:nx;
    t = (k-1)/fs;
    h = plot(t,x,'k',t,xf,'r',t,y,'b');
    set(h(1),'LineWidth',0.8);
    set(h(1),'Color',[0.8 0.8 0.8]);
    set(h(2),'LineWidth',1.5);
    set(h(3),'LineWidth',1.5);
    set(gca,'Layer','Top');
    box off;
    xlabel('Sample');
    ylabel('Signal');
    if exist('AxisSet','file'),
        AxisSet;
        end;    
    legend('Input Signal','Linear Filter Output','Nonlinear Filter Output',4);
    zoom on;
    end;
    
%====================================================================
% Process Output Arguments
%====================================================================      
if nargout==0,
    clear y;
    return;
    end;    