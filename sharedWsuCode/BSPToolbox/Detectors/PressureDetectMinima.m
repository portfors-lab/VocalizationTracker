function [f] = PressureDetectMinima(x,fsa,lfa,hfa,pfa)
%PressureDetectMinima: Detects the starts of the beats 
%
%   f = PressureDetectMinima(x,fs,lf,hf,pf);
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default = 125 Hz      
%   lf   Heart rate low frequency bound (Hz). Default = 1 Hz
%   hf   Heart rate upper frequency bound (Hz). Default = 3.5 Hz
%   pf   Plot flag: 0=none (default), 1=screen
%
%   f    Percusion peak (P1) index, samples
%
%   Heuristic DSP beat detection algorithm for pressure 
%   signals based on spectral heart rate estimation, rank 
%   order detection feature detection, derivative and 
%   inflection point estimation, nearest neighbor decicion 
%   logic, and interbeat interval misdetection/false positive
%   correction.
%
%   The algorithm employs three bandpass filters with 
%   different cutoff frequencies as a preprocessing step 
%   to the spectral heart rate estimation, derivative 
%   calculation, and rank basedpercussion peak detection. 
%   The PSD estimate in the heart rate estimation is obtained
%   using the Blackman-Tukey method combinedmwith the harmonic
%   spectrogram. Percentile based peak detection is combined 
%   with derivative estimation to determine the candidate 
%   peaks. The interbeat interval correction logic is based 
%   on rank order statistics nonlinear filters to determine 
%   the misdetections and false positives. 
%
%   Example: Find the percussion peaks in the ICP signal;
%
%      load ICP; 
%      PressureDetectMinima(icp, fs);
%
%   Version 0.00.00.00 MA
%
%   References: TBC
%
%   See also PressureDetect, PressureDetectRank, and ECGDetectQRS.

%==================================================================
% Variables
% x   = input signal
% fs  = sampling frequency
% lf  = heart rate low frequency boundary
% hf  = heart rate high frequency boundary
% Lx  = length of x
% Ldx = length of dx (derivative of x)
% wl  = window length
% xf1 = digital filtering 1
% xf2 = digital filtering 2
% xf3 = digital filtering 3
% mhr = median heart rate estimate
% f0  = feature level 0: all maxima in x
% f1  = feature level 1: rank-based feature detection in x
% f2  = feature level 2: rank-based inflection point detection in x
% f3  = featrue level 3: nearest-neighbor feature detection
%===================================================================
%x = icpcomposite;
%fs = 125;
%lf = 1;
%hf = 3.5;
%pf = 1;

%===========================================================================
% Process function arguments
%===========================================================================
if nargin<1 | nargin>5,
    help PressureDetectMinima;
    return;
end;

fs = 125;                                        % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
end;

lf = 1;                                          % Defult lower f bound, Hz
if exist('lfa') & ~isempty(lfa),
    lf = lfa;
end;

hf = 3.5;                                       % Defult upper f bound, Hz
if exist('hfa') & ~isempty(hfa),
    hf = hfa;
end;
        
pf = 0;                                          % Default - no plotting
if nargout==0,                                   % Plot if no output arguments
    pf = 1;
end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
end;

%================================================
% Process Input Arguments
%================================================
x  = x(:);
Lx = length(x);
wl = 10*fs;

%================================================
% Developer Parameters
%================================================
% 1 = Low pressure morphology
% 2 = High pressure morphology
% 3 = Unspecified morphology (default)
% t (type of morphology)
%t = 3; 
%if t == 1;
%    rthr = 80;
%elseif t == 2;
%    rthr = 50;
%elseif t == 3;
%    rthr = 60;
%end;
rthr = 15;
    
%================================================
% Process Input Arguments
%================================================
Lx = length(x);

%================================================
% Maxima Detection: detect all maxima in x
%================================================
f0 = DetectMinima(x, [], 0);

%================================================
% Digital Filtering 1 (for Heart Rate Estimation)
%================================================
xf1 = LowPass(x,fs,3*hf,1,2,0); 
xf1 = HighPass(x,fs,(1/2)*lf,2,0);

%================================================
% Spectrum-based Heart Rate Estimation
%================================================
hr = [];
for k=1:wl:Lx                               
    if ((k+wl-1 < Lx-wl))
        xs       = xf1(k:k+wl-1);
        [p,f]    = BlackmanTukey(xs, fs);
        [hp, fh] = HarmonicPSD(p,f);
        i        = find(fh>lf & fh<=hf);
        fh       = fh(i);
        hp       = hp(i);
        [m, i]   = max(hp);
        hr       = [hr; fh(i)];
     else
        xs       = xf1(k:Lx);
        [p,f]    = BlackmanTukey(xs, fs);
        [hp, fh] = HarmonicPSD(p,f);
        i        = find(fh>lf & fh<=hf);
        fh       = fh(i);
        hp       = hp(i);
        [m, i]   = max(hp);
        hr       = [hr; fh(i)];
     end;
 end;

mhr    = median(hr);          % Median Heart Rate
T_est  = (1./hr)*fs;          % Period Estimate  
T_hat  = (1/mhr)*fs;          % Median Period
T_hatl = round(T_hat/1.75);   % Low Limit 
T_hatu = round(1.75*T_hat);   % Upper Limit

%================================================
% Digital Filtering 2 (for Derivative Calculation)
%================================================
xf2  = lowpass(x,fs,3*mhr,1,2,0);
xf2  = highpass(xf2,fs,(1/2)*lf,2,0);

%================================================
% Derivative Calculation & Maximum Slope Points
%================================================
dx   = diff(xf2);             % Derivative
Ldx  = length(dx);
f2   = [];
for k=1:wl:Ldx-1,
    if ((k+wl-1 <= Ldx-wl)),
        ds     = dx(k:k+wl-1);
        thrds  = prctile(ds,90);
        f2s    = detectmaxima(ds,thrds);
        f2     = [f2; f2s+k-1];
    else
        ds     = dx(k:Ldx);
        thrds  = prctile(ds,90);
        f2s    = detectmaxima(ds,thrds);
        f2     = [f2; f2s+k-1];
    end;
end;

%================================================
% Digital Filtering 3 (for Rank-Based Feature Det)
%================================================
%xf3 = lowpass(x,fs,10*mhr,1,2,0);
xf3 = highpass(x,fs,(1/2)*lf,2,0);
f1  = [];
for k=1:wl:Lx-1,
    if ((k+wl-1 <= Lx-wl)),
        xs    = xf3(k:k+wl-1);
        thr   = prctile(xs,rthr);
        f1s   = DetectMinima(xs,thr);
        f1    = [f1;f1s+k-1];
    else
        xs    = xf3(k:Lx);
        thr   = prctile(xs,rthr);
        f1s   = DetectMinima(xs,thr);
        f1    = [f1;f1s+k-1];
    end;
end;

%================================================
% Nearest Neighbor Feature Detection
%================================================
f3  = [];
Lf2 = length(f2);
for k=1:Lf2,
    I       = find(f1<f2(k));
    [m,f3i] = min(abs(f1(I)-f2(k)));
    f3      = [f3; f1(I(f3i))];
end;
f3  = unique(f3);
f3  = sort(f3);
Lf3 = length(f3);

%f3t = [];
%for k=1:Lf3,
%    [m,f3ti] = min(abs(f0-f3(k)));
%    f3t      = [f3t; f0(f3ti)];
%end;
%f3  = f3t;


%================================================
% IBI Correction
%================================================
ibi = diff(f3);
md  = find(ibi > T_hatu);
Lmd = length(md);
ni  = f3(md);
nf  = f3(md+1);
fc  = [];
for k=1:Lmd,
    x1  = xf2(ni(k):nf(k));
    x2  = xf3(ni(k):nf(k));
    dx1 = diff(x1);
    f2c = DetectMaxima(dx1,prctile(dx1,90));
    f1c = DetectMinima(x2, prctile(x2,rthr));
    ftmp = [];
    for k2=1:length(f2c),
        i       = find(f1c < f2c(k2));
        [m, i2] = min(abs(f1c(i)-f2c(k2)));
        ftmp    = [ftmp; f1c(i(i2))];
    end;
    fc = [fc; ftmp+ni(k)];
end;
f = [f3;fc];
f = unique(f);
f = sort(f);

ibi2 = diff(f);
y10  = RankOrder(ibi2,fs,10);
y90  = RankOrder(ibi2,fs,90);
y50  = MedianFilter(ibi2,fs);

md90 = find(ibi2 > y90);
ni   = f(md90);
nf   = f(md90+1);
Lmd90= length(md90);
fc  = [];
for k=1:Lmd90-1,
    x2  = xf3(ni(k):nf(k));
    f1c = DetectMinima(x2, prctile(x2,rthr),0); 
    fc  = [fc; f1c+ni(k)];
end;
f = [f;fc];
f = unique(f);
f = sort(f);

ibi2 = diff(f);
y1   = RankOrder(ibi2,5*fs,1);
y95  = RankOrder(ibi2,5*fs,95);
y50  = MedianFilter(ibi2,fs);
md95 = find(ibi2 > y95);
ni   = f(md95);
nf   = f(md95+1);
Lmd95= length(md95);
fc   = [];
for k=1:Lmd95,
    x2  = xf3(ni(k):nf(k));
    f1c = DetectMinima(x2, prctile(x2,rthr),0);
    fc = [fc; f1c+ni(k)];
end;
f = [f;fc];
f = unique(f);
f = sort(f);

ibi3  = diff(f);
ifp   = find(ibi3 < T_hatl);
ifp   = ifp+1;
Lifp  = length(ifp);
cnt   = 0;
while ~isempty(ifp),
    fp    = f(ifp(1));
    f     = f(f~=fp);
    ibi3  = diff(f);
    ifp   = find(ibi3 < T_hatl);
    ifp   = ifp+1;
    cnt   = cnt+1;
    if cnt == 5*Lifp; break; end;
end;

ibi4 = diff(f);
y1   = RankOrder(ibi4,5*fs,1);
y99  = RankOrder(ibi4,5*fs,99);
md99 = find(ibi4 > y99);
fp1  = find(ibi4 < y1);

%================================================
% Plotting
%================================================
if pf == 1,
    t  = (1:Lx)./fs;
    tf = f./fs;
    figure;
    figureset(1, 'wide');
    h1 = axes('Position', [0.10 0.32 0.80 0.60]);     
    h  = plot(t,x,'b',tf,x(f),'r.');
    set(h1, 'Xtick', []);
    set(h, 'Markersize', 15);
    ylabel('Pressure Signal & Detection');
    axis tight;
    axisset(8);
    box off;
    
    h2   = axes('Position', [0.10 0.10 0.80 0.20]);
    f2f  = (diff(f)-1)/fs;
    f2fi = (f(1:length(f)-1)+f(2:length(f)))/2;
    h    = plot(f2fi, f2f, f2fi, f2f, 'b');
    ylabel('IBI,s');
    xlabel('Time, s');
    ylim([T_hatl/fs T_hatu/fs]);
    axis tight;
    axisset(8);
    box off;
end;

%=====================================================================
% Take care of outputs
%=====================================================================
%if nargout==0,
%    clear('f');
%end; 



        