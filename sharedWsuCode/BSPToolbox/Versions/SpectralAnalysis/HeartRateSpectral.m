function hr = HeartRateSpectral(x,fsa,wla,pfa);
%HeartRateSpectral: Estimate the heart rate (hr)
%
%   [hr] = HeartRateSpectral(x,fs,wl,pfa)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz).  Default=125 Hz      
%   wl   Window length (s). Default = 10 s, wl must be greater than 10 
%   pf   Plot fla: 0=none (default), 1=screen
%
%   hr   heart rate estimate 
%
%   Estimates the heart rate of the input signal. The input signal is 
%   expected to be a human physiologic signal such as ICP, ABP, or ECG.  
%   The output hr is column vector containing the heart rate time series. 
%   The function estimates the heart rate within each window using the PSD.
%
%   Example: Estimate the HR of an intracranial pressure signal 
%   sampled at 125 Hz and plot the results:
%
%      load ICP; 
%      [hr] = HeartRateSpectral(icp,fs,10,1);
%
%   Version  1.00 MA
%
%   See also EstimatePSD, Blackman and Periodogram.


%==================================================================
% Process function arguments
%==================================================================
if nargin<1 | nargin>4,
    help HeartRateSpectral;
    return;
    end;

fs = 125;                                 % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa), 
    fs = fsa;
    end;
 
wl = 10*fs;                               % Default window
if exist('wla') & ~isempty(wla),
    wl = wla*fs;
    end;
    
pf = 0;                                  % Default - no plotting
if nargout==0,                           % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%==================================================================
% Decimate Input
%==================================================================
xi = x(:);                                  % Original input, xi
ft = 100;                                   % New sampling frequency
R  = round(fs/ft);                          % Decimation factor
xid= decimate(x,R);                         % Decimate by  R
x  = decimate(x,R);                         % Decimate by  R
fs = fs/R;                                  % New fs for detection
LD = length(x);                             % Length of data

%=================================================================
% Remove Trend
%=================================================================
xlp = lowpass(x, fs, 0.5, 1, 2, 0);
x   = x - xlp;

%=================================================================
% Estimate PSD and Heart Rate
%=================================================================
hr       = [];
for k=1:wl:LD                               
    if ((k+wl-1 < LD-wl))
        xs       = x(k:k+wl-1);
        [p,f]    = EstimatePSD(xs, fs);
        [hp, fh] = HarmonicPSD(p,f);
        i        = find(fh>0.5 & fh<=3.5);
        fh       = fh(i);
        hp       = hp(i);
        [m, i]   = max(hp);
        hr       = [hr; fh(i)];
     else
        xs       = x(k:LD);
        [p,f]    = EstimatePSD(xs, fs);
        [hp, fh] = HarmonicPSD(p,f);
        i        = find(fh>0.5 & fh<=3.5);
        fh       = fh(i);
        hp       = hp(i);
        [m, i]   = max(hp);
        hr       = [hr; fh(i)];
     end;
 end;
 

 
%==================================================================
% Plotting 
%==================================================================
if pf==1
    figure(1)
    figureset(1)
    subplot(2,1,1)
    AxisSet(7,'Comic Sans MS');
    plot((1:length(x))./fs, xid)
    title('Input Time Series');
    ylabel('Amplitude');
    xlabel('Time, s');
    axis tight;
    xl = xlim;
    grid on;

    subplot(2,1,2)
    h = plot(hr);
    AxisSet(7,'Comic Sans MS');
    set(h, 'Linewidth', 1.5);
    title('Heart Rate Estimation');
    xlabel('Index');
    ylabel('Heart Rate, Hz');
    grid on;
end;
 
 