function [m] = SpectralHRV(p,f,pfa)
%SpectralHRV: Standard spectral heart rate variability metrics.
%
%   [m] = SpectralHRV(p,f,pf);
%
%   p    Estimated Power Spectral Density of NN intervals (NNI)
%        or instantaneous heart rate (IHR)
%   f    Frequencies at which p is estimated (Hz).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   m    Vector of the four HRV metrics recommended by the task 
%        force: [HF,LF,LFN,HFN,LHR,VLF,ULF,LLS,TP].
%
%   This function calculates the eight spectral heart rate 
%   variability (HRV) metrics recommended by the task force:
%      
%      LF    Power in LF range (0.04 - 0.15 Hz)
%      HF    Power in HF range (0.15 - 0.40 Hz)
%      LFN   Ratio of LF to (TP - VLF)
%      HFN   Ratio of HF to (TP - VLF)
%      LHR   Ratio of LF to HF
%      VLF   Power in VLF range (0.003 - 0.04 Hz)
%      ULF   Power in ULF range (<= 0.003 Hz)
%      LLS   Slope of least squares regression line of log of PSD
%            versus log of frequency
%      TP    Total power (variance)
%
%   Example: Calculate the HRV metrics of the first five minutes
%   of an ECG.
%
%      load ECG;
%      x   = ecg(1:fs*5*60);
%      di  = ECGDetectREnergy(x,fs);
%      nr  = length(di);
%      rri = diff(di)/fs;
%      t   = (di(2:nr)-1)/fs;
%      fsi = 5;
%      ti  = 0:(1/fsi):5*60; 
%      y   = SmoothSeries(t,rri,ti,0.5);
%      [p,f] = BlackmanTukey(y,fsi,5*60,5000);
%      m   = SpectralHRV(p,f,1);
%
%   Task Force of the European Society of Cardiology and the American 
%   Society of Pacing and Electrophysiology, "Heart rate variability:
%   Standards of measurement, physiological interpretation, and 
%   clinical use," Circulation, vol. 93, no. 5, pp. 1043-1065, Aug. 
%   1996.
%
%   Version 0.00.00.23 JM
%
%   See also Detectors, Metrics, and SmoothSeries.

%====================================================================
% Error Checking
%====================================================================    
if nargin<2,
    help SpectralHRV;
    return;
    end;
    
%====================================================================
% Process function arguments
%====================================================================       
pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing
%====================================================================  
p  = p(:); % Convert to column vector
f  = f(:); % Convert to column vector
np = length(p);

%====================================================================
% Calculate Metrics
%====================================================================  
vhi = find(f>0.400);
hfi = find(f>0.150 & f<=0.400);
lfi = find(f>0.040 & f<=0.150);
vli = find(f>0.003 & f<=0.040);
uli = find(f~=0    & f<=0.003);
nz  = find(f>0     & f<=0.400);

hfp = sum(p( hfi))*2/(2*np+1);
lfp = sum(p( lfi))*2/(2*np+1);
vlp = sum(p( vli))*2/(2*np+1);
ulp = sum(p( uli))*2/(2*np+1)+p(1)/(2*np+1);
tp  = sum(p(2:np))*2/(2*np+1)+p(1)/(2*np+1);
lfn = lfp/(tp-vlp);
hfn = hfp/(tp-vlp);
lhr = lfp/hfp;

fl = log(f(nz));
pl = log(p(nz));
A  = [ones(length(fl),1) fl];
w  = pinv(A)*pl;
lls = w(2);

%====================================================================
% Postprocessing
%====================================================================    
m = [lfp,hfp,lfn,hfn,lhr,vlp,ulp,lls,tp];

%====================================================================
% Plot Log-Log PSD w/Slope
%====================================================================
if pf==1,
    figure;
    FigureSet(1);
    h = plot(f(uli),p(uli),'m',f(vli),p(vli),'g',f(lfi),p(lfi),'r',f(hfi),p(hfi),'b',f(vhi),p(vhi),'k');
    set(h,'LineWidth',1.2);
    hold on;
        h = plot(0.400*[1 1],[0 1.05*max(p)],'k:');
        h = plot(0.150*[1 1],[0 1.05*max(p)],'k:');
        h = plot(0.040*[1 1],[0 1.05*max(p)],'k:');
        h = plot(0.003*[1 1],[0 1.05*max(p)],'k:');
        hold off;
    xlim([0 0.50]);
    ylim([0 1.05*max(p)]);
    title('Spectral HRV');
    xlabel('Frequency (Hz)');
    ylabel('PSD');
    box off;
    AxisSet;
    legend('ULF','VLF','LF','HF');
    
    figure;
    FigureSet(2);
    fp = [min(fl) max(fl)];
    pp = fp*w(2) + w(1);
    fp = exp(fp);
    pp = exp(pp);
    h = loglog(f(nz),p(nz),'b',fp,pp,'r');
    set(h,'LineWidth',1.2);
    xlim([min(f(nz)) max(f(nz))]);
    ylim([min(p(nz)) max(p(nz))]);
    title('Spectral HRV (Log-Log)');
    xlabel('Frequency (Hz)');
    ylabel('PSD');
    box off;
    AxisSet;
    end
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('m');
    end;
