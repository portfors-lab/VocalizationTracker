function [y] = TimeFrequencySynthesis (X, fsa, wla, pfa)
%TimeFrequencySynthesis: Time-Frequency Synthesis function
%
%   [y] = TimeFrequencySynthesis(X,fs,wl,pf)
%
%   X       Input signal(must be a matrix)
%   fs      Sampling frequency (Hz). Default=1. 
%   wl      Length of window to use (sec). Default = 1024 samples.
%           If a vector, specifies entire window.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   Y       Time domain signal
%
%   Synthesizes a frequency-domain signal (such as the spectrogram)
%   back to its time domain representation by performing the inverse
%   FFT, and dividing by the specified window. When using this 
%   function to reverse a matrix output by the Spectrogram function, 
%   several factors must be taken into consideration: 1) the sampling
%   rate input (fs) must be the same as was used when calling the 
%   Spectrogram; 2) The window length entered as an input to the
%   Spectrogram must be an integer power of two minus one (eg. 256-1
%   =255), and the corresponding window length to the TimeFrequency-
%   Synthesis function must be the next integer power of two (256);
%   3) The number of time pixels entered to the Spectrogram must be
%   equal to the length of the signal. Note that since the Spectrogram
%   subtracts the mean of the signal as one of its preprocessing steps,
%   the time-domain signal recovered with the TimeFrequencySyntesis
%   function is centered at zero.
%
%   Example: Create a spectrogram of the ICP signal and bring the signal
%   back to the time domain and plot it.
%
%      load ICP.mat;
%      icpd = decimate(icp, 15);
%      [S,t,f] = Spectrogram(icpd,125/15,blackman(255),[],[],length(icpd),1);
%      [y] = TimeFrequencySynthesis(S,125/15,256,1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, 1997.
%
%   Version 1.00 CC
%
%   See also Spectrogram.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>4,
    help TimeFrequencySynthesis;
    return;
    end;

fs = 1;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;  

wl = wla;

pf = 0; % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%--------------------------------------------------------------------
% Preprocessing
%--------------------------------------------------------------------sx = size(X);           % Size of input matrix

sx = size(X);
vl = sx(1);             % Vertical length of the input
hl = sx(2);             % Horizontal length of the input

Y  = zeros(hl,wl);
wn = blackman(wl);
w  = [wn(2);wn(2:wl-1);wn(wl-1)];
cy = 1;

%--------------------------------------------------------------------
% Time Frequency Synthesis
%--------------------------------------------------------------------
for i = 1:hl
    Z  = [X(1:vl,i);conj(X(vl-1:-1:2,i))];
    zw = real(ifft(Z));
    z  = zw./w;
    Y(hl-i+1,:) = z';
end

for j = -hl+1:1:wl-1
    y(cy) = median(diag(Y,j));
    cy = cy +1;
end

y = y(round((wl-1)/2):length(y)-round((wl-1)/2));


%--------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------
if exist('pf') & pf == 1,
    figure
    plot((0:length(y)-1)./fs, y);
    title('Time Domain Reconstruction');
    xlabel('Seconds');
    ylabel('Amplitude');
    AxisSet(10); 
    FigureSet(1);
end


