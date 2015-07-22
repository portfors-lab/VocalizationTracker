function [hp,fh] = HarmonicPSD(p,fa,nha,tha,pfa)
%HarmonicPSD: Estimate the total (harmonic) PSD versus frequency.
%
%   [hp,fh] = HarmonicPSD(p,f,th,pf)
%
%   p    Estimated power spectral density (PSD)  
%   f    Vector of frequencies at which p was estimated at. 
%        Default = 0 to almost 1.0 Hz (equivalent fs = 2 Hz).
%   nh   Number of harmonics to add. Default = 5.
%   th   Threshold for fraction of harmonic power added. Default = 2.     
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   hp   Estimated harmonic PSD (column vector).
%   fh   Frequencies at which estimate was made.
%
%   Estimates the total power at each frequency including the power
%   of the harmonics. This method was described in detail in the
%   reference below.
%
%   Example: Estimate the harmonic PSD of an ECG signal signal 
%   sampled at 500 Hz using 11 harmonics and plot the results.
%
%      load ECG;
%      [p,f] = BlackmanTukey(ecg-mean(ecg),fs,20);
%      HarmonicPSD(p,f,11);
%
%   J. McNames, C. Crespo, M. Aboy, J. Bassale, L. Jenkins, 
%   B. Goldstein, "Harmonic Spectrogram for the Analysis of 
%   Semi-periodic Physiologic Signals," submitted to 
%   EMBS and BME 2002. 
%
%   Version 1.00 JM
%
%   See also SPECTRUM, WINDOW, and SpectralAnalysis.

%====================================================================
% Process function arguments, fill in defaults
%====================================================================
if nargin<1,
    help HarmonicPSD;
    return;
    end;

np = length(p);
f  = (0:np-1)/(np); % Default
if exist('fa') & ~isempty(fa),
    f = fa;
    end;    
    
nh = 5; % Default value   
if exist('nha') & ~isempty(nha),
    nh = round(max(1,nha));
    end;

th = 2; % Default value    
if exist('tha') & ~isempty(tha),
    th = max(0,tha);
else 
    th = 2; % Default value
    end;
    
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
if nh==1,
    h = p; % If only using the fundamental, then hp = p;
    return;
    end;
    
%====================================================================
% Main Loop
%==================================================================== 
np = length(p);      % Length of input PSD
fi = 1:floor(np/nh); % Frequency indices of harmonic PSD estimate, hp
nf = length(fi);     % Number of frequencies of hp
fh = f(fi);          % Frequencies of hp
hp = zeros(nf,1);    % Allocation for hp
for c = 1:nf,
    id = (fi(c)-1)*(1:nh) + 1;
    pr              = p(id);
    pr(pr>th*pr(1)) = th*pr(1); % Apply threshold
    hp(c)           = sum(pr);    
    end;
    
%====================================================================
% Postprocessing
%====================================================================     
hp = hp(:); % Make into a column vector
fh = fh(:); % Make into a column vector

%====================================================================
% Plot Results
%====================================================================     
if pf==1,      
    figure;
    FigureSet(1);
    subplot(2,1,1);
        h = plot(f(fi),p(fi));
        set(h,'LineWidth',1.2);
        ylim([0 max(hp)*1.05]);
        xlim([0 max(f(fi))]);
        box off;
        ylabel('PSD');
        AxisSet;
    subplot(2,1,2);
        h = plot(fh,hp);
        set(h,'LineWidth',1.2);
        ylim([0 max(hp)*1.05]);
        xlim([0 max(f(fi))]);
        box off;
        xlabel('Frequency (Hz)');
    ylabel('Harmonic PSD');
    AxisSet;
    end;

    
if nargout==0,
    clear('hp','fh');
    return;
    end;
