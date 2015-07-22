function [p,f] = ModifiedPeriodogram(x,fsa,wna,nfa,pfa)
%ModifiedPeriodogram: Estimate PSD using a modified periodogram.
%
%   [p,f] = ModifiedPeriodogram(x,fs,wn,nf,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wn   Window to use. Default = Blackman. Must be same length as x.
%   nf   Number of frequencies to evaluate.
%        Default = max(128,round(length(x)/2)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   p    Power spectral density. 
%   f    Frequencies at which p is estimated (Hz).
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using a Modified Periodogram. This method estimates the PSD by 
%   calculating the FFT of the windowed signal. The window 
%   multiplication in the time domain is equivalent to convolution 
%   (a form of smoothing) in the frequency domain and trades variance 
%   of the PSD estimate for increased bias. Note that, like the 
%   periodogram, this is not a consistent estimator and the 
%   variability of the estimate cannot be controlled by any of the 
%   user-specified parameters.
%
%   Since the PSD is an estimate of power, multiplication by the 
%   window in the time domain is equivalent to convolution by the
%   square of the window in the frequency domain. The mean is removed 
%   from the signal prior to estimation to prevent an impulse at 0 Hz 
%   from dominating the smoothed estimate.
%
%   The estimated PSD is scaled such that Parseval's relation is 
%   approximately approximately satisfied: 
%                          +pi
%      var(x) ~= inv(2*pi) int p(w) dw ~= sum(p)/length(p).
%                          -pi
%
%   Example: Estimate the PSD of an electrocardiogram signal and 
%   plot the results. 
%
%      load NoisyECG.mat;
%      x  = decimate(ecg(1:50e3),5);
%      fs = fs/5;
%      ModifiedPeriodogram(x,fs,[],5000);
%
%   M. Hayes, Statistical Digital Signal Processing and Modeling. 
%   New York: John Wiley & Sons, 1996, pp. 408-412.
%
%   J. G. Proakis, C. M. Rader, F. Ling, C. L. Nikias, M. Moonen, 
%   and I. K. Proudler, Algorithms for Statistical Signal Processing.
%   Saddle River, NJ: Prentice Hall, 2002, pp. 447-448.
%
%   Version 1.00 JM
%
%   See also SPECTRUM, WINDOW, and SpectralAnalysis.


%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help ModifiedPeriodogram;
    return;
    end;

nx = length(x);
if nx==0,
    error('Signal is empty.');
    end;
    
%====================================================================
% Process function arguments
%====================================================================
fs = 1;                              % Default sampling rate
if exist('fsa') & ~isempty(fsa)
    fs = fsa;
    end

wn = blackman(nx);                   % Default window shape
if exist('wna') & ~isempty(wna) & length(wna)==nx,
    wn = wna(:); 
    end;
        
nf = max(128,nx);                    % Default No. of frequencies 
if exist('nfa') & ~isempty(nfa)
    nf = nfa;
    end   
    
pf = 0;                              % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Preprocessing and Memory Allocation
%====================================================================
x   = x (:).';               % Convert to a row vector
wn  = wn(:).';               % Convert to a row vector
wn  = wn/sqrt(mean(wn.^2));  % Normalize window to make asymptotically unbiased 
mx  = mean(x);               % Signal mean
vx  = var(x);                % Signal variance
x   = x  - mean(x);          % Remove the signal mean

%====================================================================
% Initialize Variables
%====================================================================
nz   = 2^(ceil(log2(nf))+1); % No. of zeros to pad with
fi   = 1:(nz/2+1);           % Index of frequencies to evaluate at
f    = fs*((fi-1)/nz)';      % Frequencies
nf   = length(f);            % No. of frequencies
fr   = fs/nz;                % Frequency resolution (frequency difference between adjacent estimates)

%====================================================================
% Estimate PSD
%====================================================================
xf  = fft(x.*wn,nz);    % Calculate FFT of windowed signal   
xf  = xf(fi);           % Truncate negative frequencies
p   = abs(xf).^2;       % Calculate PSD estimate
p   = p/nx;             % Scale PSD estimate

%====================================================================
% Postprocessing
%====================================================================  
f = f(:);               % Ensure is column vector
p = p(:);               % Ensure is column vector

%====================================================================
% Plot PSD
%====================================================================
if pf==1,
    figure;
    FigureSet;
    plot(f,p);
    title('Modified Periodogram Estimated PSD');
    xlabel('Frequency (Hz)');
    ylabel('PSD (scaled)');
    xlim([0 fs/2]);
    ylim([0 1.02*max(p)]);
    box off;
    zoom on;
    AxisSet;
    end

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('p','f');
    end;

    

