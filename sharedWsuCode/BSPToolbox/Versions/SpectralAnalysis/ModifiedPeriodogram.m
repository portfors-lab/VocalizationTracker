function [psd,f] = ModifiedPeriodogram(x,fsa,wla,tpa,pfa)
%ModifiedPeriodogram: Estimate the PSD using the averaging modified 
%                     periodogram
%
%   [psd,f] = ModifiedPeriodogram(x,fs,w,tp,pf)        
%
%   x     Input signal       
%   fs    Sampling frequency in Hz. Default=125      
%   wl    Length of window to use (sec). If a vector, specifies 
%         entire window 
%         Default = 1001 samples
%   tp    Portion of serie to be tapered. Default=50
%   pf    Plot format: 0=none, 1=screen. Default=0
%
%   psd   Power spectral density 
%   f     Frequencies (Hz) at which p is estimated 
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the modified periodogram method. The PSD of the input 
%   signal is calculated by taking the absolute value of the FFT 
%   of the windowed signal squared. The user can specify the type 
%   of window (blackman, bartlett, or hanning), the sampling  
%   frequency, and whether or not to plot the results. In addition, 
%   the user can specify what portion of the serie should be tapered  
%   to reduce leakage. 
%
%   If w is a scalar, then it specifies the window type as shown 
%   above. If w is a vector, then w is treated as a user-specified 
%   window. If it is a different length than x, zero padding is used 
%   to make the windowing possible.
%
%   Example: Estimate the PSD of an intracranial blood pressure 
%   signal sampled at 125 Hz using a Blackman window and plot the 
%   results
%
%      load ICP.mat; 
%      ModifiedPeriodogram(icp, 125, 2);
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.408-412, 1996.
%
%   Version 1.00 JB
%
%   See also BlackmanTukey, EigenVectorPSD, HeartRate,
%   MaximumEntropyPSD, MinimumNormPSD, Music, Welch, and Wigner.


%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>5,
    help ModifiedPeriodogram;
    return;
    end;
    
if ( size(x,1)== 1)
    x = x';
end


fs = 125;                          % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;

wl = 1001;          % Default window length
wn = blackman(wl)'; % Default window shape
if exist('wla') & ~isempty(wla),
    if length(wla)==1, % User specified window length only
        wl = max(3,round(wla*fs));
        if ~rem(wl,2),
            wl = wl + 1; % Make odd
            end;
        wn = blackman(wl)'; % Default window shape
    else % User specified the actual window, not just the length
        wn = wla(:)';
        wl = length(wn);
        if ~rem(wl,2),
            wn = [wn 0]; % Append a zero to make odd
            wl = wl + 1;
            end;
        end;
    end;

tp  = 0.5;                           % Default tapering
if exist('tpa') & ~isempty(tpa),
    tp = tpa;
    end;

pf  = 0;                             % Default - no plotting
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end; 
   
 
  
% ---------------------------------------------------------------------
% Detrending and de-meaning.
% ---------------------------------------------------------------------
    nwx          = x(:) - mean(x(:));
    n            = length(nwx);

% ---------------------------------------------------------------------
% Taper data
% ---------------------------------------------------------------------
    m = floor(n * tp);
    t = (1:m)-0.5;
    t1 = (n-m+1:n)+0.5;
    W(1:m) = .5*[1-cos(pi*t./m)];
    W(m+1:n-m) = 1;
    W(n-m+1:n) = .5*[1-cos(pi*t1./m)];
    x= x.*W';
 
    
% ---------------------------------------------------------------------
% Periodogram Parameters
% ---------------------------------------------------------------------
    n  = max(n, length(wl));
    xx           = zeros(1,n);
    xx(1:length(x)) = x;
    x               = xx;

    n2  = ceil(log2(n));
    n3 = 0:n-1;
    N  = 2^n2;
  


% ---------------------------------------------------------------------
% Estimate PSD
% ---------------------------------------------------------------------
    xw    = x.*wl'/norm(wl);
    PX    = abs(fft(x, N)).^2/N;
    hr    = n3*(1/(N*(1/fs)));
    f     = hr(1:N/2);
    psd   = PX(1:N/2);
    psd   = psd';
    f     = f';
% ---------------------------------------------------------------------
% Plotting (TBC by James)
% ---------------------------------------------------------------------
if (pf ==1 | nargout == 0)
    figure;
    FigureSet(1);
   stem(f, abs(psd));grid
   title('Modified Periodogram ');
   xlabel('Hz');
   ylabel('PSD');
end
