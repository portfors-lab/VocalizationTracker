function [S,t,f,handles] = NonparametricSpectrogram(x,fs,varargin)
%NonparametricSpectrogram: Sliding-window periodogram of signal
%
%   [S,t,f, handles] = NonparametricSpectrogram(x,fs,windowLength,
%                           frequencyRange,nFrequencies,
%                           nTimes,plotType);
%   [S,t,f] = NonparametricSpectrogram(x,fs,'Parameter',Value,...);
%
%   x               Input signal.
%   fs              Sample rate (Hz).
%   windowDuration  Length of data window to use (sec). Default = 
%                   1024 samples. If a vector, specifies window itself.
%   autocorrelationDuration
%                   Duration of autocorrelation window to use (sec).
%                   Default = infinite, in which case the Spectrogram
%                   uses the modified periodogram.
%   frequencyRange  Minimum and maximum frequencies to display (Hz).
%                   Default = [0 fs/2].
%   nFrequencies    Number of frequencies to evaluate (vertical 
%                   resolution). Default = max(128,round(nWindow/2)).
%   nTimes          Requested number of times (horizontal pixels) 
%                   to evaluate. Default = min(400,length(x)).
%   plotType        Plot type: 0=none (default), 1=screen, 2=current 
%                   figure, 'None', Screen', 'CurrentFigure'
%   spectrogramType Specifies the format of the plotted and returned
%                   spectrogram.
%                   0 = square root of the power spectral density (default)
%                   1 = power spectral density
%                   2 = abs(power spectral density) on a decibel scale
%   dbReference     The reference power for the decibel scale.
%                   Only used when spectrogramType=2.
%                   Default = max(max(abs(power spectral density)))
%   dbRange         The dynamic range, in decibels, of the spectrogram
%                   relative to the peak calculated decibel value.
%                   Only used when spectrogramType=2. Default = 100dB.
%   timeOffset      Start time (s) of recording.
%   colorLimit      Upper color limit as percentile of the spectrogram.  
%                   Default = 97.
%   colormap        colormap to use for spectrogram
%   figure          A handles structure from previous spectrogram figure,
%                   will delete old spectrogram before plotting new one
%   visible         visiblity of figure, default='on'
%
%   S       Matrix containing the image of the Spectrogram.
%   t       Times at which the spectrogram was evaluated (s).
%   f       Frequencies at which the spectrogram was evaluated (Hz).
%   handles Structure containing the handles of the figure and components:
%           figure, haSpectrogram, haSignal.  This is the same as the input
%           parameter 'figure'
%
%   Calculates estimates of the spectral content at the specified 
%   times using a modified periodogram. The mean of the signal is 
%   removed as a preprocessing step. The square root of the power 
%   spectral density is calculated and displayed. To limit 
%   computation, decimate the signal if necessary to make the upper 
%   frequency range approximately equal to half the sampling 
%   frequency. 
%   
%   If only the window length is specified, the blackman window is 
%   used. The specified window length should be odd. If called with
%   a window with an even number of elements, a zero is appended to
%   make the window odd.
%
%   Example: Generate the spectrogram of an intracranial pressure
%   signal using a Blackman-Harris window that is 45 s in duration.
%
%      load ICP.mat; 
%      icpd = decimate(icp,15);
%      nWindow   = round(45*fs/15);
%      NonparametricSpectrogram(icpd,fs/15,blackmanharris(nWindow));
%
%   M. Hayes, "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, 1996.
%
%   Version 1.00 JM
%
%   See also spectrogram, window, decimate, Beatogram, and 
%   Cohereogram.

% Modified Amy Boyle Aug 2010; extracted handles from graphics
%====================================================================
% Error Checking
%====================================================================    
if nargin<2
    help NonparametricSpectrogram;
    return;
end
    
if isempty(x)
    error('Signal is empty.');
end

if 0  && var(x)==0
    error('Signal is constant.');
end
      
%====================================================================
% Calculate Basic Signal Statistics
%====================================================================    
nx    = length(x);
xRaw  = x;
xMean = mean(x);
    
%====================================================================
% Process Function Arguments
%====================================================================
nWindow    = min(1024,nx);                                 % Default window length
window     = blackman(nWindow)';                           % Default window shape
nAutocorrelation = inf;                              % Autocorrelation window duration
windowAutocorrelation  = nan;
fmin       = 0;                                            % Lowest frequency to display
fmax       = fs/2;                                         % Highest frequency to display
nf         = max(256,round(nWindow/2));                    % Number of frequencies
nSteps     = min(400,nx);                                  % Number of times
pf         = 0;                                            % No plot by default
specType   = 0;                                            % sqrt(PSD) by default
dbRange    = 100;                                          % Default dynamic range of dB scales
timeOffset = 0;
colorLimitPercentile = 97;
visibility = 'on';
cmap = 'jet';
if nargout==0, pf = 1; end;                                % Plot if no output arguments

if(nargin > 2 && ischar(varargin{1}))
    if mod(nargin-2,2) ~= 0
        error(['Unexpected number of arguments. Optional arguments ' ...
              'must either all be specified as property value pairs or ' ...
              'all be specified in the order listed in the documentation.']);
    end
    for c1 = 1:2:nargin-2
        if(~ischar(varargin{c1}))
            error(['Error parsing arguments: Expected property name string at argument ' num2str(c1+1)]);
        end
        value = varargin{c1+1};
        switch(lower(varargin{c1}))
            case {'windowlength','windowduration'},
                if length(value)==1,                       % User specified window length only
                    nWindow = max(3,round(value*fs));
                    if ~rem(nWindow,2), nWindow = nWindow + 1; end; % Make odd
                    window = blackman(nWindow)';           % Default window shape
                else                                       % User specified the actual window, not just the length
                    window = value(:)';                    % Make into a row vector
                    nWindow = length(window);
                    if ~rem(nWindow,2),
                        window = [window 0];                       % Append a zero to make odd
                        nWindow = nWindow + 1;
                    end
                end
            case {lower('AutocorrelationDuration')},
                if length(value)==1                        % Specified duration only
                    nAutocorrelation = max(3,round(value*fs));
                    if ~rem(nAutocorrelation,2), nAutocorrelation = nAutocorrelation + 1; end; % Make odd
                    autocorrelationWindow = blackman(nAutocorrelation).'; 
                else                                       % User specified the actual window, not just the length
                    autocorrelationWindow  = value(:)';
                    nAutocorrelation = length(autocorrelationWindow);
                    if ~rem(nAutocorrelation,2),
                        autocorrelationWindow = [autocorrelationWindow 0]; % Append a zero to make odd
                    end
                    nAutocorrelation = length(autocorrelationWindow);
                end
                
            case 'frequencyrange'
                if length(varargin{c1+1})==1,
                    error('Frequency range must be a 2-element vector.');
                end;
                fmin = max(varargin{c1+1}(1),0);
                fmax = min(varargin{c1+1}(2),fs/2);
            case 'nfrequencies'
                nf = varargin{c1+1};
            case 'ntimes'
                nSteps = min(varargin{c1+1},nx);
            case 'plottype'
                pf = varargin{c1+1};
            case 'spectrogramtype'
                specType = varargin{c1+1};
            case 'dbreference'
                dbRef = varargin{c1+1};
            case 'dbrange'
                dbRange = varargin{c1+1};
            case lower('timeOffset')
                timeOffset = varargin{c1+1};                
            case lower('ColorLimit')
                colorLimitPercentile = value;          
            case 'figure'
                handles = varargin{c1+1};
            case 'visible'
                visibility = varargin{c1+1};
            case 'colormap'
                cmap = varargin{c1+1};
            otherwise
                error(['Unrecognized property: ''' varargin{c1} '''']);
        end
    end
else
    warning('Use newer style of function arguments with (ParameterName, ParameterValue) pairs. Old style will be removed in future versions.');
    switch(nargin-2)
        case 0,
        case 1,
            [wla] = varargin{:};
        case 2,
            [wla fra] = varargin{:};
        case 3,
            [wla fra nfa] = varargin{:};
        case 4,
            [wla fra nfa nsa] = varargin{:};
        case 5,
            [wla fra nfa nsa pfa] = varargin{:};
        otherwise
            error('Too many arguments');
    end

    if exist('wla','var') && ~isempty(wla),
        if length(wla)==1,                                     % User specified window length only
            nWindow = max(3,round(wla*fs));
            if ~rem(nWindow,2),
                nWindow = nWindow + 1;                                   % Make odd
            end;
            window = blackman(nWindow)';                            % Default window shape
        else                                                         % User specified the actual window, not just the length
            window = wla(:)';
            nWindow = length(window);
            if ~rem(nWindow,2),
                window = [window 0];                                   % Append a zero to make odd
                nWindow = nWindow + 1;
            end;
        end;
    end;

    if exist('fra','var') && length(fra)==1,
        error('Frequency range must be a 2-element vector.');
    end;
    if exist('fra','var') && ~isempty(fra),
        fmin = max(fra(1),0);
        fmax = min(fra(2),fs/2);
    end;        

    if exist('nfa','var') && ~isempty(nfa),
        nf = nfa;
    end;   

    if exist('nsa','var') && ~isempty(nsa),
        nSteps = min(nsa,nx);
    end;   
    if exist('pfa','var') && ~isempty(pfa),
        pf = pfa;
    end;
end;

    
%====================================================================
% Preprocessing
%====================================================================    
window = window(:).';                                      % Make into a row vector
x      = x (:).';                                          % Make into a row vector
x      = x - xMean;                                        % Remove sample average

%====================================================================
% Variable Allocations & Initialization
%====================================================================
nWindow        = length(window);
nSamplesStep   = (nx-1)/(nSteps-1);                        % Step size (samples)
nSamplesStep   = max(nSamplesStep,1);
iTime          = 1:nSamplesStep:nx;                        % Evaluation times (in samples)
nSteps         = length(iTime);                            % No. evaluation times
nSamplesPadded = (nf-1)*(fs/(fmax-fmin));                  % No. window points needed in FFT
nSamplesPadded = 2^(ceil(log2(nSamplesPadded)));           % Convert to power of 2 for FFT
if isinf(nAutocorrelation) 
    bBlackmanTukey = false;
else                              
    bBlackmanTukey         = true; 
    nAutocorrelation = min([nAutocorrelation;nWindow]);
    if nAutocorrelation > nSamplesPadded/2-1
        warning('Autocorrelation window length (%d) too large for FFT length (%d). Reduced to maximum size.',nAutocorrelation,nSamplesPadded); %#ok
        nAutocorrelation = nSamplesPadded/2-1;
    end
end
if ~bBlackmanTukey && nWindow>nSamplesPadded 
    warning('Length of FFT (%d) smaller than segment (%d). Using Blackman Tukey estimate to compensate.',nSamplesPadded,nWindow); 
    nAutocorrelation = nSamplesPadded/2-1;
    bBlackmanTukey   = true;
end

b0   = floor(nSamplesPadded*(fmin/fs))+1;                 % Lower frequency bin index 
b1   = ceil (nSamplesPadded*(fmax/fs))+1;                 % Upper frequency bin index
fi   = b0:b1;                                              % Frequency indices
f    = fs*(fi-1)/nSamplesPadded;                           % Frequencies
nf   = length(fi);                                         % No. of frequencies that PSD is evaluated at 

%====================================================================
% Main loop
%====================================================================
S  = zeros(nf,nSteps);                                     % Power Spectral Density
for c1 = 1:nSteps,
    ic  = iTime(c1);                                       % Sample index of segment center
    i0  = round(ic-nWindow/2);                             % Index of first segment element               
    i1  = i0 + (nWindow-1);      % Index of last segment element
    %k   = (i0:i1);                                        % Collection of indices
    nZerosPrepend    = max(-i0+1,0);                       % Number of zeros to append at front
    nZerosAppend     = max(i1-nx,0);                       % Number of zeros to append at end
    %k   = k(1+nZerosPrepend:nWindow-nZerosAppend);        % Subset of valid indice    
    k                = max(i0,1):min(i1,nx);
    xSegment         = [x(1)*ones(1,nZerosPrepend) x(k) x(nx)*ones(1,nZerosAppend)];                
	xSegmentWindowed = xSegment.*window;
    if bBlackmanTukey
        nSamplesPaddedSegment = 2^nextpow2(2*nWindow-1); % Convert to power of 2 for FFT
        xSegmentDft           = fft(xSegmentWindowed,nSamplesPaddedSegment);                                            % Calculate FFT of x
        xAutocorrelation      = ifft(xSegmentDft.*conj(xSegmentDft)); % Calculate biased autocorrelation
        xAutocorrelation      = real(xAutocorrelation(1:nWindow)).'; % Eliminate superfluous zeros and nearly zero imaginary part
        xAutocorrelation      = [xAutocorrelation(nAutocorrelation:-1:2);xAutocorrelation(1:nAutocorrelation)]; % Make two-sided
        windowAutocorrelation = blackman(length(xAutocorrelation)+2);
        windowAutocorrelation = windowAutocorrelation(2:end-1);
        xAutocorrelation      = xAutocorrelation.*windowAutocorrelation;        
        xDft                  = sqrt(abs(fft(xAutocorrelation,nSamplesPadded))); % Window autocorrelation and calculate the estimate
    else
        xDft = fft(xSegmentWindowed,nSamplesPadded);       % Discrete fourier transform
    end
    S(:,c1) = xDft(fi).';
    drawnow;
end;
t = (iTime-0.5)/fs;

%====================================================================
% Postprocessing
%====================================================================  
x = xRaw;   
t = t(:);                                                  % Convert to column vector
f = f(:);                                                  % Convert to column vector    
if specType == 1
    %Convert to power spectral density
    S = abs(S).^2;
elseif specType == 2
    %Convert to dB scale
    if ~exist('dbRef','var')
        dbRef = max(max(abs(S)));
    end
    %Turn off warnings for log(0) calculations
    warning off
    Sdb = 20*log10(abs(S)/dbRef);
    warning on
    dbMax = max(max(Sdb));
    Sdb = max(Sdb,dbMax-dbRange);
end
        

%====================================================================
% Plot Results
%====================================================================  
if pf>=1,
    if pf~=2,
        if exist('handles','var')
            delete([handles.haSpectrogram, handles.haSignal, ...
                handles.haColorbar, handles.haPowerSpectralDensity]);
        else
            handles.figure = figure('numbertitle', 'off', 'Visible', visibility);
            %FigureSet(1);
        end;
%         colormap(jet(256));
        colormap(cmap);
    end;
    
    td = 1;                                                % Time Divider
    if max(t)>2000,
        td = 60;                                           % Use minutes
    end;
    
    Tmax = (nx-0.5)/fs;
        
    fmax = max(f);
    k    = 1:length(x);
    tx   = (k-0.5)/fs;    
    tp   = t;
    tmax = max(t);
    tw   = [(nWindow/2-1)/fs,(tmax-(nWindow/2)/fs)];
    
    %--------------------------------------------------------------------------
    % Determine the Units of the Axis Labels
    %--------------------------------------------------------------------------
    xLabel   = 'Time (s)';                                     % X-axis label
    if max(timeOffset+t)>2000,                                        % Convert to minutes
        tx         = tx/60;
        tp         = tp/60;
        tmax       = tmax/60;
        tw         = tw/60;
        timeOffset = timeOffset/60;
        xLabel         = 'Time (min)';
    end    
    
    yLabel = 'Frequency (Hz)';
    if max(f)>=40e3
        f      = f/1000;
        fmin   = fmin/1000;
        fmax   = fmax/1000;
        yLabel = 'Frequency (kHz)';
    end

    %--------------------------------------------------------------------------
    % Spectrogram
    %--------------------------------------------------------------------------
    haSpectrogram = axes('Position',[0.15 0.21 0.81 0.69]);
    handles.haSpectrogram = haSpectrogram;
    if specType == 0 || specType == 1
        s = reshape(abs(S),nf*nSteps,1);  
        %p = [0 max(s)];
        p = [0 prctile(s,colorLimitPercentile)];
        imagesc(timeOffset+tp,f,abs(S),p);
    else
        imagesc(timeOffset+tp,f,Sdb);
    end

    xlim(timeOffset+[0 tmax]);
    ylim([fmin fmax]);
    %set(haSpectrogram,'YTickLabel',[]);    
    set(haSpectrogram,'YAxisLocation','Right');
    set(haSpectrogram,'XAxisLocation','Top');
    set(haSpectrogram,'YDir','normal');
    hold on;
        h = plot([1;1]*(tw+timeOffset),[0 fs],'k');
        set(h,'LineWidth',2.0);
        h = plot([1;1]*(tw+timeOffset),[0 fs],'w');
        set(h,'LineWidth',1.0);
    hold off;

    %--------------------------------------------------------------------------
    % Colorbar
    %--------------------------------------------------------------------------
    %haColorbar = axes('Position',[0.135 0.21 0.01 0.69]);    
    haColorbar = colorbar('Position',[0.135 0.21 0.01 0.69]);    
    handles.haColorbar = haColorbar;
    %colorbar(haColorbar);    
    set(haColorbar,'Box','Off')
    set(haColorbar,'YTick',[]);

    %--------------------------------------------------------------------------
    % Power Spectral Density
    %--------------------------------------------------------------------------
    haPowerSpectralDensity = axes('Position',[0.08 0.21 0.05 0.69]);
    handles.haPowerSpectralDensity = haPowerSpectralDensity;
    if specType == 0
        psd = mean((abs(S).^2).').^(1/2);
        h = plot(psd,f,'r'); 
        xlim([0 max(psd)]);        
    elseif specType == 1
        psd = mean((abs(S)).');
        h = plot(psd,f,'r'); 
        xlim([0 max(psd)]);  
    else
        psd = mean((abs(S).^2).').^(1/2);
        warning off
        psd = 20*log10(psd/dbRef);
        warning on
        psd = max(psd,dbMax-dbRange);
        h = plot(psd,f,'r'); 
        xlim([min(psd) max(psd)]);
    end
    set(h,'LineWidth',1.0);
    ylabel(yLabel);
    set(haPowerSpectralDensity,'XTick',[]);
    ylim([fmin fmax]);

    %--------------------------------------------------------------------------
    % Signal
    %--------------------------------------------------------------------------
    handles.haSignal = axes('Position',[0.15 0.10 0.81 0.10]);

    sts  = max(round(nSamplesStep/25),1);
    k    = 1:sts:nx;
    h    = plot(timeOffset+tx(k),x(k));
    set(h,'LineWidth',1.0);
    ymin = min(x);
    ymax = max(x);
    yrng = ymax-ymin;
    ymin = ymin - 0.02*yrng;
    ymax = ymax + 0.02*yrng;
    xlim(timeOffset + [0 tmax]);
    ylim([ymin ymax]);
    xlabel(xLabel);
    ylabel('Signal');
    
    %--------------------------------------------------------------------------
    % Final Settings
    %--------------------------------------------------------------------------
    %axes(haSpectrogram);
    set(gcf, 'currentaxes', haSpectrogram);  
    AxisSet(8);
    drawnow;
    frequencyLink = linkprop([haSpectrogram,haPowerSpectralDensity],'ylim');
    timeLink = linkprop([haSpectrogram,handles.haSignal],'xlim');
    setappdata(haPowerSpectralDensity,'graphics_linkprop',frequencyLink); 
    setappdata(handles.haSignal,'graphics_linkprop',timeLink); 
end
 
%====================================================================
% Process Return Arguments
%====================================================================
if specType == 2
    S = Sdb;
    clear Sdb
end
if nargout==0,
    clear('S','t','f');
end