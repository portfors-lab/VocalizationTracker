function [D,t,f] = Delayogram(x,y,fsa,wla,cla,fra,nfa,nsa,hfa,pfa);
%Delayogram: Nonstationary estimate the group delay versus time
%
%   [D,t,f] = Delayogram(x,y,fs,sl,wl,sn,ol,fr,nf,ns,pf);
%   
%   x    Input signal.
%   y    Output signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1024 samples.
%        If a vector, specifies entire window.
%   cl   Length of correlation window to use in Blackman-Tukey
%        spectral estimates. Default = 2*wl/5. 
%   fr   Minimum and maximum frequencies to display (Hz).
%        Default = [0 fs/2].
%   nf   Number of frequencies to evaluate (vertical resolution). 
%        Default = max(128,round(wl/2)).
%   ns   Requested number of times (horizontal pixels) to evaluate 
%        Default = min(400,length(x)).
%   hf   Histogram equalization: 0=none (default), 1=uniform
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   D    Matrix containing the image of the group delay vs. time.
%   t    Times at which the cohereogram was evaluated (s).
%   f    Frequencies at which the cohereogram was evaluated (Hz).
%
%   This function estimates the group delay of a transfer function
%   and plots the result as an image. 
%
%   Zero padding is used to improve frequency resolution. The 
%   group delay is calculated by using a first-order difference of
%   the estimated cross-power spectral density.
%
%        C = angle(Pxy).
%
%   If no output arguments are specified, the image is generated with 
%   a plot of the signals in the bottom axis. The color map is scaled
%   from 0 to the 99th percentile. The signal is extrapolated at the
%   edges by repeating the value of the edge.
%
%   Example:  Plot the group delay of 20 minute segments of ABP and ICP
%   data. Limit the frequency axis to only show values between 0 and 
%   4 Hz. Use a segment length of 1 minute and 10 second subsegments.
%   Use a triangular window.
%
%      load ABPICP.mat
%      x1 = decimate(abp,15);
%      x2 = decimate(icp,15);
%      fs = fs/15;
%      Delayogram(x1,x2,fs,60);
%
%   Version 1.00 JM
%
%   See also visualization.

%====================================================================
% Error checking
%====================================================================
if nargin<2,
    help Delayogram;
    return;
    end;

nx = min(length(x),length(y));
if nx==0,
    fprintf('Error: Input signals must have finite length.\n');
    return;
    end;
   
if var(x)==0,
    fprintf('Error: Signal x is constant.\n');
    return;
    end;
    
if var(y)==0,
    fprintf('Error: Signal y is constant.\n');
    return;
    end; 
   
%====================================================================
% Process function arguments & Fill in defaults
%====================================================================
fs = 1;                     % Default: 1 Hz
if exist('fsa') & ~isempty(fsa),
    fs = abs(fsa);
    end;   

wl = min(1024,nx);          % Default window length
if exist('wla') & ~isempty(wla),
    wl = max(3,round(wla*fs));
    end;
    
cl = max(3,ceil(wl/5));    % Default window length
if exist('cla') & ~isempty(cla),
    cl = min(round(cla*fs),wl);
    cl = cl + (1-rem(cl,2)); % Make cl odd
    end;
    
fmin = 0;                   % Default: 0 Hz 
fmax = fs/2;                % Default: fs/2 Hz
if exist('fr') & ~isempty(fr) & length(fr)==2,
    fmin = max(fr(1),0);
    fmax = min(fr(2),fs/2);
    end;
    
nf = max(200,round(wl/2));           % Default No. vertical pixels
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
    end;   
    
ns = min(200,nx);           % Default: 400 horizontal pixels
if exist('nsa') & ~isempty(nsa),
    ns = min(nsa,nx);
    end;         
        
hf = 0;                                 % Default - no equalization
if exist('hfa') & ~isempty(hfa),
    hf = hfa;
    end;    
    
pf = 0;                     % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing
%====================================================================    
xr = x;
yr = y;

mx = mean(x);
my = mean(y);

x  = x(:).' - mx; % Make into a row vector and remove mean
y  = y(:).' - my; % Make into a row vector and remove mean

%====================================================================
% Variable Allocations & Initialization
%====================================================================
st   = (nx-1)/(ns-1);         % Step size (samples)
st   = max(st,1);
te   = 1:st:nx;               % Evaluation times (in samples)
nt   = length(te);            % No. evaluation times
nz   = nf*(fs/(fmax-fmin));   % No. window points needed in FFT
nz   = max(2*wl+1,nz);        % Ensure cross-correlation is not aliased
nz   = 2^(ceil(log2(nz)));    % Convert to power of 2 for FFT
b0   = floor(nz*(fmin/fs))+1; % Lower frequency bin index 
b1   = ceil (nz*(fmax/fs))+1; % Upper frequency bin index
fi   = b0:b1;                 % Frequency indices
f    = fs*(fi-1)/nz;          % Frequencies
nf   = length(fi);            % No. of frequencies that PSD is evaluated at 

%====================================================================
% Main loop
%====================================================================
D  = zeros(nf,nt);            % Power Spectral Density
for c1 = 1:nt,
    ic  = te(c1);             % Sample index of segment center
    i0  = round(ic-wl/2);     % Index of first segment element               
    i1  = i0 + (wl-1);        % Index of last segment element
    k   = (i0:i1);            % Collection of indices
    a0  = max(-i0+1,0);       % Number of zeros to append at front
    a1  = max(i1-nx,0);       % Number of zeros to append at end
    k   = k(1+a0:wl-a1);      % Subset of valid indice    
    
    xp  = [x(1)*ones(1,a0) x(k) x(nx)*ones(1,a1)];        
    yp  = [y(1)*ones(1,a0) y(k) y(nx)*ones(1,a1)];        
            
	sx  = fft(xp,nz);          % Discrete fourier transform
	sy  = fft(yp,nz);          % Discrete fourier transform
    sxy = sx.*conj(sy);
    axy = unwrap(angle(sxy));    
    
    crw  = real(ifft(sxy).');
    cwd  = zeros(size(crw));
    
    hw   = (cl-1)/2;
    wd   = blackman(cl);
    cwd(1:hw+1        ) = crw(1:hw+1        ).*wd(hw + (1:hw+1));
    cwd(nz-(hw-1:-1:0)) = crw(nz-(hw-1:-1:0)).*wd(1:hw);
    
    ssxy = fft(cwd,nz);
    saxy = unwrap(angle(ssxy));
        
    dxy = diff(axy)*(nz/(2*pi));
    dxy = [dxy(1) dxy];
    
    sdxy = diff(saxy)*(nz/(2*pi));
    sdxy = [sdxy(1);sdxy];
    dxy  = sdxy.';
    dxy  = dxy/fs; % Convert from samples to seconds 
    
    if 0,    
        FigureSet(2);
        subplot(3,1,1);
            h = plot(f,axy(fi),'r',f,saxy(fi),'b');
            set(h,'Marker','.');
            ylim([-pi pi]);
            ylabel('Angle');
        subplot(3,1,2);
            h = plot(f,dxy(fi),'r',f,sdxy(fi),'b');
            set(h,'Marker','.');       
            ylim([-20 20]);
            ylabel('Delay');
        subplot(3,1,3);
            t = (-wl:wl)/fs;
            cxy   = [crw(nz-(wl-1:-1:0));crw(1:wl+1)];
            scxy  = [cwd(nz-(wl-1:-1:0));cwd(1:wl+1)];
            plot(t,cxy,'r',t,scxy,'b');            
            ylabel('Cross-Correlation');
        fprintf('Mean Raw:%5.2f Smoothed:%5.2f\n',mean(dxy(fi)),mean(sdxy(fi)));
        fprintf('Pausing...\n');
        pause;
        end;
    
    D(:,c1) = dxy(fi).';    
	end;
t = (te-1)/fs;

%fprintf('Mean overall delay: %5.3f\n',mean(mean(D)));

%====================================================================
% Postprocessing
%====================================================================  
x = xr; 
y = yr;
t = t(:); % Convert to column vector
f = f(:); % Convert to column vector    

%====================================================================
% Plot Results
%====================================================================  
if pf==1,
    td = 1; % Time Divider
    if max(t)>2000,
        td = 60; % Use minutes
        end;
    
    Tmax = (nx-1)/fs;
    
    figure;
    FigureSet(1);

    fmax = max(f);
    k    = 1:length(x);
    tx   = (k-1)/fs;    
    tp   = t;
    tmax = max(t);
    tw   = [(wl/2-1)/fs,(tmax-(wl/2)/fs)];
    xl   = 'Time (s)'; % X-axis label
    if max(t)>2000, % Convert to minutes
        tx   = tx/60;
        tp   = tp/60;
        tmax = tmax/60;
        tw   = tw/60;
        xl   = 'Time (min)';
        end    
    
    %----------------------------------------------------------------
    % Delayogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.10 0.21 0.75 0.69]);
    d = reshape(D,nf*nt,1);
    
    if hf==1, % Apply histogram equalization
        nc = 256;
        colormap(jet(nc));
        [ds,id] = sort(d);
        ns = length(d);
        i0 = 1;
        for c1 = 1:nc,
            i1 = ceil(c1*ns/nc);
            ds(i0:i1) = c1;
            i0 = i1;
            end;
        d(id) = ds;
        D = reshape(d,nf,nt);
        image(tp,f,D);
    else
        p = [prctile(d,02) prctile(d,98)];
        imagesc(tp,f,D,p); %,[Smin Smax]);
        end;
    xlim([0 tmax]);
    ylim([fmin fmax]);
    set(ha1,'YTickLabel',[]);    
    set(ha1,'XAxisLocation','Top');
    set(ha1,'YDir','normal');
    hold on;
        h = plot([1;1]*tw,[0 fs],'k');
        set(h,'LineWidth',2.0);
        h = plot([1;1]*tw,[0 fs],'w');
        set(h,'LineWidth',1.0);
        hold off;     
    
    %----------------------------------------------------------------
    % Colorbar
    %----------------------------------------------------------------    
    ha2 = colorbar;    
    set(ha2,'Box','Off')
    set(ha2,'YTick',[]);
    set(ha2,'Position',[0.97 0.10 0.02 0.80]);
    set(ha1,'Position',[0.14 0.21 0.82 0.69]);
    
    %----------------------------------------------------------------
    % Power Spectral Density
    %----------------------------------------------------------------
    ha3 = axes('Position',[0.08 0.21 0.05 0.69]);
    md  = mean(D.');
    plot(md,f); %,[Smin Smax]);
    ylim([fmin fmax]);
    xlim([min(md) max(md)]);
    ylabel('Frequency (Hz)');
    set(gca,'XAxisLocation','Top');    
    
    %----------------------------------------------------------------
    % Output Signal
    %----------------------------------------------------------------      
    ha4 = axes('Position',[0.14 0.15 0.82 0.05]);
    h    = plot(tx,y);
    set(h,'LineWidth',1.5);
    ymin = min(y);
    ymax = max(y);
    yrng = ymax-ymin;
    ymin = ymin - 0.02*yrng;
    ymax = ymax + 0.02*yrng;
    xlim([0 tmax]);
    ylim([ymin ymax]);
    ylabel('Output');
    
    %----------------------------------------------------------------
    % Input Signal
    %----------------------------------------------------------------      
    ha4 = axes('Position',[0.14 0.10 0.82 0.05]);

    h    = plot(tx,x);
    set(h,'LineWidth',1.5);
    ymin = min(x);
    ymax = max(x);
    yrng = ymax-ymin;
    ymin = ymin - 0.02*yrng;
    ymax = ymax + 0.02*yrng;
    xlim([0 tmax]);
    ylim([ymin ymax]);
    xlabel(xl);
    ylabel('Input');
    
    axes(ha1);
    AxisSet(8);
    end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('D','t','f');
    end;
