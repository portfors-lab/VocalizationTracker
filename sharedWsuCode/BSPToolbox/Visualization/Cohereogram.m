function [C,t,f] = Cohereogram(x1,x2,fsa,sla,wla,sna,ola,fra,nfa,nsa,pfa);
%Cohereogram: Nonstationary estimate the coherency versus time
%
%   [C,t,f] = Cohereogram(x1,x2,fs,sl,wl,sn,ol,fr,nf,ns,pf);
%   
%   x1   First input signal.
%   x2   Second input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   sl   Length of signal segments used to generate estimate (sec). 
%        Default = (signal duration)/100.
%   wl   Length of each of the windows applied to subsets of each 
%        segment (s). Default: sl/20. If vector, wl is used as
%        the window. Otherwise a Blackman window is applied.
%   sn   Signal to noise ratio used to bias the coherence. 
%        Default = inf.
%   ol   Overlap of windows applied to each subset (%). Default = 50%.
%   fr   Minimum and maximum frequencies to display (Hz).
%        Default = [0 fs/2].
%   nf   Number of frequencies to evaluate (vertical resolution). 
%        Default = max(128,round(wl/2)).
%   ns   Requested number of times (horizontal pixels) to evaluate 
%        Default = min(400,length(x)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   C    Matrix containing the image of the cohereogram.
%   t    Times at which the cohereogram was evaluated (s).
%   f    Frequencies at which the cohereogram was evaluated (Hz).
%
%   This function estimates the coherence (not the coherency 
%   spectrum) of two input signals in a moving window and plots the 
%   result as an image.  The coherence is a measure of correlation of 
%   estimated spectral content of two signals as a function of 
%   frequency. The range of values is 0 to 1. 
%
%   Zero padding is used to improve frequency resolution. The 
%   coherence is calculated by estimating the following ratio.
%
%        C = abs(Pxy)/sqrt(Pxx.*Pyy).
%
%   If no output arguments are specified, the image is generated with 
%   a plot of the signals in the bottom axis. The color map has a 
%   fixed range of 0 to 1. The signal is extrapolated at the edges by 
%   repeating the value of the edge.
%
%   Example:  Plot the coherence of 20 minute segments of ABP and ICP
%   data. Limit the frequency axis to only show values between 0 and 
%   4 Hz. Use a segment length of 1 minute and 10 second subsegments.
%   Use a triangular window.
%
%      load ABPICP.mat
%      x1 = decimate(abp,15);
%      x2 = decimate(icp,15);
%      fs = fs/15;
%      wn = triang(round(10*fs));
%      Cohereogram(x1,x2,fs,60,wn);
%
%   Version 1.00 JM
%
%   See also Coherency, NonparametricSpectrogram, and cohere.

%====================================================================
% Error checking
%====================================================================
if nargin<2,
    help Cohereogram;
    return;
    end;

nx = min(length(x1),length(x2));
if nx==0,
    fprintf('Error: Input signals must have finite length.\n');
    return;
    end;
   
if var(x1)==0,
    fprintf('Error: Signal x1 is constant.\n');
    return;
    end;
    
if var(x2)==0,
    fprintf('Error: Signal x2 is constant.\n');
    return;
    end; 
   
%====================================================================
% Process function arguments & Fill in defaults
%====================================================================
fs = 1;                     % Default: 1 Hz
if exist('fsa') & ~isempty(fsa),
    fs = abs(fsa);
    end;   

sl = round(nx/100);         % Default: 1000 samples     
sl = max(100,sl);           % Do not let fall below 100
if exist('sla') & ~isempty(sla)
    sl = round(sla*fs);
    end;
sl = min(nx,sl);            % Do not let exceed nx

wl = max(ceil(sl/20),2);    % Default: sl/20, must be at least 2
wn = blackman(wl);          % Default window type
if exist('wla') & ~isempty(wla),
    if length(wla)==1,      % User-specified window length
        wl = round(wla*fs); % Convert from seconds to samples
        wl = min(wl,round(sl/2)); % Cannot exceed sl/2
        wl = max(wl,2);     % Must be at least 2
        wn = blackman(wl);  % Default window shape
    else                    % User specified window
        wl = length(wla);
        wn = wla;
        end;
    end;

sn = inf;                   % Default: infinite signal to noise ration (SNR)
if exist('sna') & ~isempty(sna),
    sn = sna;                                              % User-specified (must be a positive number)
    end;
    
ol = 0.5;                   % Default: 50% overlap
if exist('ola') & ~isempty(ola),
    ol = ola/100;           % User-specified (convert from percentage to fraction)
    ol = min(ol,1);         % Must not exceed 100%
    ol = max(ol,0);         % Must be non-negative
    end;
    
fmin = 0;                   % Default: 0 Hz 
fmax = fs/2;                % Default: fs/2 Hz
if exist('fra') & ~isempty(fra) & length(fra)==2,
    fmin = max(fra(1),0);
    fmax = min(fra(2),fs/2);
    end;
    
nf = max(200,round(wl/2));           % Default No. vertical pixels
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
    end;   
    
ns = min(200,nx);           % Default: 400 horizontal pixels
if exist('nsa') & ~isempty(nsa),
    ns = min(nsa,nx);
    end;         
        
pf = 0;                     % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Variable Allocations & Initialization
%====================================================================
bs   = wl - ol*wl;            % Subsegment step size (fractional samples)
bs   = max(1,bs);             % Must be at least 1 sample 
bn   = 1+floor((sl-wl)./bs);  % No. of subsegment steps to take
bn   = max(2,bn);             % Take at least 2 steps
st   = (nx-1)/(ns-1);         % Segment step size (fractional samples)
st   = max(st,1);             % Must not be less than 1 sample
te   = 1:st:nx;               % Segment evaluation times (in samples)
nt   = length(te);            % No. evaluation times
nz   = nf*(fs/(fmax-fmin));   % No. window points needed for specified frequency resolution
nz   = 2^(ceil(log2(nz)));   % Convert to power of 2 for FFT
b0   = floor(nz*(fmin/fs))+1; % Lower frequency bin index 
b1   = ceil (nz*(fmax/fs))+1; % Upper frequency bin index
fi   = b0:b1;                 % Frequency indices
f    = fs*(fi-1)/nz;          % Frequencies
nf   = length(fi);            % No. of frequencies that PSD is evaluated at 

%====================================================================
% Preprocessing and Memory Allocation
%====================================================================
x1  = x1(:).'; % Convert to a row vector
x2  = x2(:).'; % Convert to a row vector    
wn  = wn(:).'; % Convert to a row vector

wn  = wn*sqrt(wl*inv(nz*sum(wn.^2))); % Ensure window preserves energy 

m1  = mean(x1);
m2  = mean(x2);

x1  = x1 - mean(x1);
x2  = x2 - mean(x2);

C   = zeros(nf,nt); 
S12 = zeros(ns,nf);  % Estimated cross power of x1 and x2
S11 = zeros(ns,nf);  % Estimated power of each subsegment of x1
S22 = zeros(ns,nf);  % Estiamted power of each subsegment of x2

sw  = blackman(ns+2);
sw  = sw(2:ns+1);
sw  = sw*ns/sum(sw);
SW  = sw*ones(1,nf); % Segment window for weighted correlation estimate

%====================================================================
% Main Loop
%====================================================================
for c1 = 1:nt, % Loop over segments  
    ic  = te(c1);                   % Sample index of segment center
    i0  = round(ic-sl/2);           % Index of first segment element               
    i1  = i0 + (sl-1);              % Index of last segment element
    k   = (i0:i1);                  % Collection of indices
    a0  = max(-i0+1,0);             % Number of zeros to append at front
    a1  = max(i1-nx,0);             % Number of zeros to append at end
    k   = k(1+a0:sl-a1);            % Subset of valid indices
            
    s1 = [x1(1)*ones(1,a0), x1(k), x1(nx)*ones(1,a1)]; % Segment of x1
    s2 = [x2(1)*ones(1,a0), x2(k), x2(nx)*ones(1,a1)]; % Segment of x2
                
    for c2 = 1:bn, % Loop over subwindows
        i0  = round((c2-1)*bs + 1); % Initial segment index for subsegment
        i1  = i0 + (wl-1);          % Ending segment index for subsegment
        i1  = min(i1,sl);           % Must not exceed segment length
        i0  = i1-(wl-1);            % Adjust initial index, if necessary    
        k   = (i0:i1);
                    
        %disp(round([i0 i1 sl]));
        X1  = fft(s1(k).*wn,nz);
        X2  = fft(s2(k).*wn,nz);
         
        S12(c2,:) = X1(fi).*conj(X2(fi));
        S11(c2,:) = abs(X1(fi)).^2;
        S22(c2,:) = abs(X2(fi)).^2;    
        
        end;        
        
    S12 = S12.*SW; % Weight the estimates for better temporal resolution
    S11 = S11.*SW;
    S22 = S22.*SW;
        
    num     = abs(sum(S12)).^2;
    d1      = sum(S11);
    d2      = sum(S22);
    if ~isinf(sn) && var(s1)>0, 
        d1 = d1+(1/sn)*sum(d1)/var(s1); % Add equivalent signal noise
    end; 
    if ~isinf(sn) && var(s2)>0, 
        d2 = d2+(1/sn)*sum(d2)/var(s2); % Add equivalent signal noise
    end; 
    den     = d1.*d2;
    id      = find(den~=0);       % Avoid divide by zero
    c       = zeros(nf,1);        % When den==0, set coherency to 0
    c(id)   = num(id)./den(id);   % Calculate coherence
    c       = sqrt(c);            % Convert coherence to coherency
    C(:,c1) = c;
    end;

%====================================================================
% Postprocessing
%====================================================================  
x1 = x1 + m1;
x2 = x2 + m2;

t = (te(:)-1)/fs; % Convert sample times (in samples) to column vector (in seconds)
f  = f(:);        % Make into a column vector   

%====================================================================
% Plot Results
%====================================================================  
if pf,
    if pf~=2,
        figure;
        FigureSet(1);
        colormap(jet(256));
        end;    
    
    fmax = max(f);
    k    = 1:nx;
    tx   = (k-1)/fs;    
    tp   = t;
    tmax = max(t);
    ts   = [(sl/2)/fs,(tmax-(sl/2)/fs)];
    xl   = 'Time (s)'; % X-axis label
    if max(t)>2000, % Convert to minutes
        tx   = tx/60;
        tp   = tp/60;
        tmax = tmax/60;
        ts   = ts/60;
        xl   = 'Time (min)';
        end    
     
    %----------------------------------------------------------------
    % Coherogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.15 0.32 0.78 0.58]);
    imagesc(tp,f,C,[0 1]); %,[Smin Smax]);
    xlim([0 tmax]);
    ylim([0 fmax]);
    set(ha1,'XAxisLocation','Top');
    set(ha1,'YDir','normal');
    set(ha1,'YTickLabel',[]);
    hold on;
        h = plot([1;1]*ts,[0 fs],'k',[1;1]*ts,[0 fs],'w');
        set(h(1:2),'LineWidth',2.0);
        set(h(3:4),'LineWidth',1.0);
        hold off;  
    title('Coherogram');
    
    %----------------------------------------------------------------
    % Colorbar
    %----------------------------------------------------------------    
    ha2 = colorbar;    
    set(ha2,'Box','Off')
    set(ha2,'Position',[0.94 0.10 0.02 0.80]);
    set(ha1,'Position',[0.14 0.32 0.79 0.58]);
    
    %----------------------------------------------------------------
    % Coherence
    %----------------------------------------------------------------
    ha3 = axes('Position',[0.08 0.32 0.05 0.58]);
    plot(mean(C'),f); %,[Smin Smax]);
    ylim([0 fmax]);   
    xlim([0 1]);
    ylabel('Frequency (Hz)');
    set(gca,'XAxisLocation','Top');
    
    %----------------------------------------------------------------
    % Signal 1
    %---------------------------------------------------------------- 
    ha4  = axes('Position',[0.14 0.21 0.79 0.10]);
    h    = plot(tx,x1);
    set(h,'LineWidth',1.5);
    ymin = min(x1);
    ymax = max(x1);
    yrng = ymax-ymin;
    ymin = ymin - 0.02*yrng;
    ymax = ymax + 0.02*yrng;
    xlim([0 tmax]);
    ylim([ymin ymax]);
    set(ha4,'xTickLabel',[]);
    ylabel('x1');
            
    %----------------------------------------------------------------
    % Signal 2
    %----------------------------------------------------------------      
    ha5 = axes('Position',[0.14 0.10 0.79 0.10]);
    h    = plot(tx,x2);
    set(h,'LineWidth',1.5);
    ymin = min(x2);
    ymax = max(x2);
    yrng = ymax-ymin;
    ymin = ymin - 0.02*yrng;
    ymax = ymax + 0.02*yrng;
    xlim([0 tmax]);
    ylim([ymin ymax]);
    ylabel('x2');
    xlabel(xl);    
    
    axes(ha1);
    AxisSet(8);
    end
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear C;
    clear t;
    clear f;
    end;