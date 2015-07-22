function [f] = PressureDetector(x,fsa,lfa,hfa,pfa)
%PressureDetector: Pressure detector algorithm (ICP, ABP)
%
%   f = PressureDetector(x,fs,lf,hf,pf);
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
%   Algorithm PRESSUREDETECTV2(x, [lf hf])
%   begin
%     F0 = Detect all maxima in x; 
%     Select a window: wl
%     for i=1 to length of x; incrementing by wl, do
%       XF1 = Bandpass-Filter x (remove f < lf; f > 3*hf);  
%       HR  = Estimate heart rate in XF1 based on HarmonicPSD
%       XF2 = Bandpass-Filter x (remove f < lf; f > 2*median(HR);
%       F2  = Find maxima in Derivative(x) using rank filter;
%       XF3 = Bandpass-Filter x (remove f < lf; f > 10*median(HR);
%       F1  = Find percussion peaks in XF3 using rank filter & maxima;
%     end
%     Nearest neighbor feature detection:
%     Among the points selected as candidate percussion peaks based 
%     on their amplitude (F1), find those which are immediately 
%     preceded by a large slope (F2)
%     for i=1 to length of x, do
%       I   = Find all F1 > F2(i);
%       F3  = Find F1 closest to F2(i);
%     end;
%     Perform IBI based correction
%     IBI = Find the interbeat interbasl in F3;
%     MD  = Find indexes where IBI is < 1/2*HR or IBI is > 2*HR
%     F   = Correct misdetection and flase positives (input: MD, F3, F0)
%  end
%
%   Example: Find the percussion peaks in the ICP signal;
%
%      load ICP; 
%      PressureDetector(icp, fs);
%
%   Version 0.00.00.00 MA
%
%   References: TBC
%
%   See also PressureDetect, PressureDetectRank, and ECGDetectQRS

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
%x = icp;
%fs = 125;
%lf = 1;
%hf = 3.5;
%pf = 1;

%===========================================================================
% Process function arguments
%===========================================================================
if nargin<1 | nargin>5,
    help PressureDetector;
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

%===========================================================================
% Process function arguments
%===========================================================================
Lx = length(x);
wl = fs*60*5;
f  = [];
for k=1:wl:Lx                              
    if ((k+wl-1 < Lx-wl))
        xs       = x(k:k+wl-1);
        fss = pressuredetectv2(xs,fs,lf,hf,0);
        f = [f;fss+k-1];
    else
        xs       = x(k-wl:Lx);
        fss = pressuredetectv2(xs,fs,lf,hf,0);
        f = [f;fss+k-wl-1];
     end;
 end;
f = unique(f);

 %================================================
% Plotting
%================================================
if pf == 1,
    t  = (1:length(x))./fs;
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
    axis tight;
    axisset(8);
    box off;
end;

%=====================================================================
% Take care of outputs
%=====================================================================
if nargout==0,
    clear('f');
end; 