function [Hf,Ha,Ht] = TrackHarmonicsNPS(x,fs,varargin)
%TrackHarmonicsNPS: Track the harmonics of a signal using a sliding-window periodogram
%
%   [Hf,Ha,Ht] = TrackHarmonicsNPS(x,fs,windowLength,
%                                  frequencyRange,nFrequencies,
%                                  nTimes,plotType);
%   [Hf,Ha,Ht] = TrackHarmonicsNPS(x,fs,'Parameter',Value,...);
%
%   x               Input signal.
%   fs              Sample rate (Hz).
%   windowDuration  Length of window to use (sec). Default = 1024 
%                   samples. If a vector, specifies entire window.
%   frequencyRange  Minimum and maximum frequencies to search for
%                   harmonic peaks in.
%                   Default = [0 fs/2].
%   nHarmonics      Number of harmonics to detect
%                   Default = 1.
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
%
%   S    Matrix containing the image of the Spectrogram.
%   t    Times at which the spectrogram was evaluated (s).
%   f    Frequencies at which the spectrogram was evaluated (Hz).
%
%   Calculates estimates of the spectral content at the specified 
%   times using a modified periodogram. The mean of the signal is 
%   removed as a preprocessing step. The square root of the power 
%   spectral density is calculated. Each periodogram calculation
%   is lowpass filtered to smooth the contour before spectral
%   peak calculation. The algorithm does not attemp to "thread"
%   together the peaks, but instead returns the frequencies and
%   amplitudes of the top nHarmonic peaks. To limit computation,  
%   decimate the signal if necessary to make the upper 
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
%      wl   = round(45*fs/15);
%      NonparametricSpectrogram(icpd,fs/15,blackmanharris(wl));
%
%   M. Hayes, "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, 1996.
%
%   Version 1.00 LH Adapted from NonperiodicSpectrogram
%
%   See also spectrogram, window, decimate, Beatogram, and 
%   Cohereogram.

%====================================================================
% Error Checking
%====================================================================    
if nargin<2
    help TrackHarmonicsNPS;
    return;
end
    
if isempty(x)
    error('Signal is empty.');
end

if var(x)==0
    error('Signal is constant.');
end
      
%====================================================================
% Calculate Basic Signal Statistics
%====================================================================    
nx = length(x);
xr = x;
mx = mean(x);
    
%====================================================================
% Process Function Arguments
%====================================================================
%Set default values
wl = min(1024,nx);                                         % Default window length
wn = blackman(wl)';                                        % Default window shape
fmin = 0;                                                  % Lowest frequency to display
fmax = fs/2;                                               % Highest frequency to display
nh = 3;                                                    % Number of harmonics
nf = max(256,round(wl/2));                                 % Number of frequencies
ns = min(400,nx);                                          % Number of times
pf = 0;                                                    % No plot by default
specType = 0;                                              % sqrt(PSD) by default
dbRange = 100;                                             % Default dynamic range of dB scales
timeOffset = 0;
frequencyCutoffRatio = 40;
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
end;

if(nargin > 2 && ischar(varargin{1}))
    if mod(nargin-2,2) ~= 0
        error(['Unexpected number of arguments. Optional arguments ' ...
              'must either all be specified as property value pairs or ' ...
              'all be specified in the order listed in the documentation.']);
    end
    for c1 = 1:2:nargin-2
        if(~ischar(varargin{c1}))
            error(['Error parsing arguments: Expected property ' ...
                   'name string at argument ' num2str(c1+1)]);
        end
        switch(lower(varargin{c1}))
            case {'windowlength','windowduration'},
                if length(varargin{c1+1})==1,              % User specified window length only
                    wl = max(3,round(varargin{c1+1}*fs));
                    if ~rem(wl,2),
                        wl = wl + 1;                       % Make odd
                    end;
                    wn = blackman(wl)';                    % Default window shape
                else                                       % User specified the actual window, not just the length
                    wn = varargin{c1+1}(:)';
                    wl = length(wn);
                    if ~rem(wl,2),
                        wn = [wn 0];                       % Append a zero to make odd
                        wl = wl + 1;
                    end;
                end;
            case 'frequencyrange'
                if length(varargin{c1+1})==1,
                    error('Frequency range must be a 2-element vector.');
                end;
                fmin = max(varargin{c1+1}(1),0);
                fmax = min(varargin{c1+1}(2),fs/2);
            case 'nharmonics'
                nh = varargin{c1+1};
            case 'nfrequencies'
                nf = varargin{c1+1};
            case 'ntimes'
                ns = min(varargin{c1+1},nx);
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
            otherwise
                error(['Unrecognized property: ''' varargin{c1} '''']);
        end
    end
else
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
            wl = max(3,round(wla*fs));
            if ~rem(wl,2),
                wl = wl + 1;                                   % Make odd
            end;
            wn = blackman(wl)';                            % Default window shape
        else                                                         % User specified the actual window, not just the length
            wn = wla(:)';
            wl = length(wn);
            if ~rem(wl,2),
                wn = [wn 0];                                   % Append a zero to make odd
                wl = wl + 1;
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
        ns = min(nsa,nx);
    end;   
    if exist('pfa','var') && ~isempty(pfa),
        pf = pfa;
    end;
end;

    
%====================================================================
% Preprocessing
%====================================================================    
wn = wn(:).';                                              % Make into a row vector
x  = x (:).';                                              % Make into a row vector
x  = x - mx;                                               % Remove sample average

%====================================================================
% Variable Allocations & Initialization
%====================================================================
st   = (nx-1)/(ns-1);                                      % Step size (samples)
st   = max(st,1);
te   = 1:st:nx;                                            % Evaluation times (in samples)
nt   = length(te);                                         % No. evaluation times
nz   = (nf-1)*(fs/(fmax-fmin));                            % No. window points needed in FFT
nz   = 2^(ceil(log2(nz)));                                 % Convert to power of 2 for FFT
b0   = floor(nz*(fmin/fs))+1;                              % Lower frequency bin index 
b1   = ceil (nz*(fmax/fs))+1;                              % Upper frequency bin index
fi   = b0:b1;                                              % Frequency indices
f    = fs*(fi-1)/nz;                                       % Frequencies
nf   = length(fi);                                         % No. of frequencies that PSD is evaluated at 

%====================================================================
% Main loop
%====================================================================
Ha  = nan(nh,nt);                                         % Harmonic Amplitudes
Hf  = nan(nh,nt);                                         % Harmonic Frequencies
for c1 = 1:nt,
    ic  = te(c1);                                          % Sample index of segment center
    i0  = round(ic-wl/2);                                  % Index of first segment element               
    i1  = i0 + (wl-1);                                     % Index of last segment element
    k   = (i0:i1);                                         % Collection of indices
    a0  = max(-i0+1,0);                                    % Number of zeros to append at front
    a1  = max(i1-nx,0);                                    % Number of zeros to append at end
    k   = k(1+a0:wl-a1);                                   % Subset of valid indice    
    y   = [x(1)*ones(1,a0) x(k) x(nx)*ones(1,a1)];        
        
	s   = y.*wn;
	dft = fft(s,nz);                                       % Discrete fourier transform
    
    periodogram = dft(fi).';

    if c1 == 2000
        foo = 1;
    end
    [peak_amplitudes peak_frequencies] = findHarmonicPeaks(abs(periodogram), ...
                                                           fs, ...
                                                           nh, ...
                                                           f, ...
                                                           fi, ...
                                                           frequencyCutoffRatio);
    Ha(:,c1) = peak_amplitudes;
    Hf(:,c1) = peak_frequencies;

end;
t = (te-0.5)/fs;

%====================================================================
% Postprocessing
%====================================================================  
x = xr;   
t = t(:);                                                  % Convert to column vector
f = f(:);                                                  % Convert to column vector    
if specType == 1
    %Convert to power spectral density
    Ha = abs(Ha).^2;
elseif specType == 2
    %Convert to dB scale
    if ~exist('dbRef','var')
        dbRef = max(max(abs(Ha)));
    end
    %Turn off warnings for log(0) calculations
    warning off
    HaDb = 20*log10(abs(Ha)/dbRef);
    warning on
    dbMax = max(max(HaDb));
    HaDb = max(HaDb,dbMax-dbRange);
end
        

%====================================================================
% Plot Results
%====================================================================  
if pf>=1,
    if pf~=2,
        figure;
        FigureSet(1);
        colormap(jet(256));
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
    tw   = [(wl/2-1)/fs,(tmax-(wl/2)/fs)];
    xl   = 'Time (s)';                                     % X-axis label
    if max(timeOffset+t)>2000,                             % Convert to minutes
        tx         = tx/60;
        tp         = tp/60;
        tmax       = tmax/60;
        tw         = tw/60;
        timeOffset = timeOffset/60;
        xl         = 'Time (min)';
    end    
    
    %----------------------------------------------------------------
    % Spectrogram
    %----------------------------------------------------------------
    haSpectrogram = axes('Position',[0.08 0.21 0.85 0.69]);
    
    plot(timeOffset+tp,Hf);

    xlim(timeOffset+[0 tmax]);
    ylim([fmin fmax]);
    ylabel('Frequency (Hz)');
    set(haSpectrogram,'XTickLabel',[]);   
    hold on;
        h = plot([1;1]*(tw+timeOffset),[0 fs],'k');
        set(h,'LineWidth',2.0);
        h = plot([1;1]*(tw+timeOffset),[0 fs],'w');
        set(h,'LineWidth',1.0);
    hold off;     
    
    %----------------------------------------------------------------
    % Signal
    %----------------------------------------------------------------      
    haSignal = axes('Position',[0.08 0.10 0.85 0.10]);

    sts  = max(round(st/25),1);
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
    xlabel(xl);
    ylabel('Signal');
        
    %----------------------------------------------------------------
    % Final Settings
    %----------------------------------------------------------------      
    axes(haSpectrogram);
    AxisSet(8);
    drawnow;
    linkaxes([haSpectrogram,haSignal],'x');
end
    
%====================================================================
% Process Return Arguments
%====================================================================
Ht = t;
if specType == 2
    Ha = HaDb;
    clear HaDb
end
if nargout==0,
    clear('Ha','Hf','t','f');
end

function [peak_amplitudes peak_frequencies] = findHarmonicPeaks(periodogram, ...
                                                                sampleRate, ...
                                                                nHarmonics, ...
                                                                freqs, ...
                                                                freqIdx, ...
                                                                frequencyCutoffRatio)

harmonicPeakDetectionNoiseFloor = 0.0;
cutoff = sampleRate/frequencyCutoffRatio;

if var(periodogram) == 0
    max_idx = [];
else
%     smoothedPeriodogram = Lowpass(periodogram,sampleRate,cutoff,4);
    smoothedPeriodogram = periodogram;
    max_idx = DetectMaxima(smoothedPeriodogram,harmonicPeakDetectionNoiseFloor,0);
end
num_peaks = length(max_idx);
max_vals = periodogram(max_idx);

peaks = zeros(nHarmonics,2);
if num_peaks > 0
    peaks2sort = [max_vals(:) max_idx(:)];
    %Sort peaks by amplitude of peaks, in decending order
    sorted_peaks = sortrows(peaks2sort,-1);
    %Prune if there are too many found
    sorted_peaks = sorted_peaks(1:min(num_peaks,nHarmonics),:);
    %Sort peaks by index
    sorted_peaks = sortrows(sorted_peaks,2);
    peaks(1:min(num_peaks,nHarmonics),:) = sorted_peaks;
end
%Add to signal length matrix to form one signal per hermonic
peak_amplitudes = nan(nHarmonics,1);
peak_frequencies = nan(nHarmonics,1);
for hIdx = 1:nHarmonics
    if peaks(hIdx,2) > 0
        peak_amplitudes(hIdx,1) = peaks(hIdx,1);
        peak_frequencies(hIdx,1) = freqs(peaks(hIdx,2));
    end
end



