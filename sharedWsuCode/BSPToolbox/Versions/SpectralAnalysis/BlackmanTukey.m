function [psd,f] = BlackmanTukey(x,fsa,wla,ma,tpa,pfa)
%BlackmanTukey: Estimate the PSD using the Blackman-Tukey method
%
%   [psd,f] = BlackmanTukey(x,fs,wl,m,tp,pf)        
%
%   x     Input signal    
%   fs    Sampling frequency in Hz, Default=125      
%   wl    Length of window to use (sec). Default=1001 
%         If a vector, specifies entire window
%   m     Length of lag window. Default=32
%   tp    Portion of serie to be tapered. Default=50
%   pf    Plot format: 0=none, 1=screen. Default=0
%
%   psd   Power spectral density 
%   f     Frequencies (Hz) at which p is estimated 
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the Blackman-Tukey method. The PSD of the input 
%   signal is estimated by convolving the signal with the FFT 
%   of the autocorrelation window whose length is determined by m. 
%   The user can specify the type of window (blackman, bartlett, 
%   or hanning), the sampling frequency, and whether or not 
%   to plot the results. 
%
%   If w is a scalar, then it specifies the window type as shown 
%   above. If w is a vector, then w is treated as a user-specified 
%   window. If it is a different length than x, zero padding is used 
%   to make the windowing possible.
%
%   Example: Estimate the PSD of an intracranial pressure signal 
%   sampled at 125 Hz using a Blackman window and plot the results.
%
%      load ICP.mat;
%      BlackmanTukey(icp, fs, blackman(1024));
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.408-412, 1996.
%
%   Version 1.00 JB
%
%   See also EigenVectorPSD, HeartRate,
%   MaximumEntropyPSD, MinimumNormPSD, ModifiedPeriodogram
%   Music, Welch, and Wigner.


% ==================================================
% Default parameters
% =================================================
if ( nargin < 1 | nargin > 6)
    help BlackmanTukey;
    break;
end;

fs = 125;                        % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;
    
m  = 32;
if exist('ma') & ~isempty(ma),
    m = ma;
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


tp  = 0.5;
if exist('tpa') & ~isempty(tpa),
    tp = tpa;
    end;

pf  = 0;                                % Default - no plotting
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end; 

    
% =================================================
% Estimate PSD
% ================================================
    n = length(x);
    n            = max(n, size(wl,1));
    xx           = zeros(1,n);
    xx(1:length(x)) = x;
    x               = xx;

    n2  = ceil(log2(n));
    N   = 2^n2;
    c = zeros(m,1);
    n = floor((n+1)/2);
    f = (0:N-1)'/(N *(1/fs));
    f = f(1:N/2);
 
    R = SetCovarMatrix(x, m);
    r = [fliplr(R(1,2:m)), R(1,1), R(1,2:m)];
    m = 2* m-1;

   
    r = r'.* wl;
    Px = abs(fft(r, N));
    Px(1) = Px(2);
    psd = Px(1:N/2);
    freq = f;
   

% ===============================================
% Plotting (TBC by James)
% ===============================================
if (pf ==1 | nargout == 0)
    figure;
    FigureSet(1);
   stem(f, abs(psd));grid
   title('Blackman-Tukey based PSD');
   xlabel('Hz');
   ylabel('PSD');
end
