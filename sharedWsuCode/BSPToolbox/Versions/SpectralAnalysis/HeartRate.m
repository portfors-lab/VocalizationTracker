function [hr] = HeartRate(x,fra,fsa,nha,wla,wl2a)
%HeartRate: Estimate the heart rate of a pressure signal
%
%   [hr] = HeartRate(x,fr,fs,nh,wl,wl2);
%
%   x     Input signal
%   fr    Minimum and maximum frequencies to display. Default=[0 fs/2]
%   fs    Sample rate (Hz). Default:1 
%   nh    Number of harmonics to include. Default:5
%   wl    Length of window to use (samples). Default=2^12
%   wl2   Specifies the number of FFT points used to calculate 
%         the PSD estimate. Default=2^15
%
%   hr    Heart Rate estimate
%
%   Calculates an estimate of the heart rate given a
%   pressure wave signal; The algorithm uses multiple harmonics
%   present in the estimate of the power spectral density of 
%   the signal to evaluate the heart rate.
%
%   Example: Generate the heart rate of an intracranial pressure
%   signal.
%      load ICP.mat; 
%      HeartRate(icp,[1 3],fs,5,2^12, 2^15);
%
%   Version 1.00 JB
%
%   See also BlackmanTukey, EigenVectorPSD
%   MaximumEntropyPSD, MinimumNormPSD, ModifiedPeriodogram
%   Music, Welch, and Wigner.

%====================================================================
% Process function arguments, fill in defaults
%====================================================================
if ( nargin < 2 | nargin > 6)
    help HeartRate;
    break;
end;

n = length(x);
if n==0,
    fprintf('ERROR: Signal is empty.\n');
    return;ds
end;

fs = 1;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
end;    

nh = 5;
if exist('nha') & ~isempty(nha),
    nh = nha; 
end;    
    
Fmin = 0;
Fmax = (fs/2)/nh;
if exist('fra') & ~isempty(fra),
    Fmin = max(fra(1),0);
    Fmax = min(fra(2),(fs/2)/nh);
end

lgth = 10;
if exist('lgtha') & ~isempty(lgtha),
    lgth = lgtha;
end


wl = 2^12;
if exist('wla') & ~isempty(wla),
    wl = wla; 
end;    

wl2 = 2^15;
if exist('wl2a') & ~isempty(wl2a),
    wl2 = wl2a; 
end;    

%====================================================================
% Estimate the heart rate
%====================================================================
R = floor(.8*fs/(2*nh*Fmax));      % Decimation coefficient 
dx = decimate(x, R);                % Decimate the signal
[psd, freq]      = pwelch(dx - mean(dx) ,blackman(wl),[],wl2,.8*fs/R);
Ind = find( (freq <= Fmax) & ( freq >= Fmin ) );
pwr = ones(length(Ind), nh);

for t = 1:length(Ind)
    for t1 = 1:nh
        e(t1) = psd(t1*Ind(t));
    end
    E(t) = sum(e);
end
[iv, i] = max(E);
hr= freq(Ind(i));


