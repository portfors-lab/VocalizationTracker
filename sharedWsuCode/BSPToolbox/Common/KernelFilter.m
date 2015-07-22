function [y,h,s] = KernelFilter(x,fs,fc,foa,pfa);
%function [y,h,s] = KernelFilter(x,fs,fc,foa,rpa,lba,mda,pfa);
%KernelFilter: Non-negative FIR filter based on kernels.
%
%   [y,h,s] = KernelFilter(x,fs,fc,fo,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz).   
%   fc   Cutoff frequency (Hz).
%   fo   Filter order (samples). Must be even. Default = 200.
%   pf   Plot flag: 0=none (default), 1=screen
%
%   y    Filtered Signal.
%   h    Filter impulse response.
%   s    Structure of filter statistics.
%
%   Creates an optimal FIR filter with non-negative impulse response
%   that is as close as possible to an ideal lowpass filter. The 
%   current implementation uses an Epanechnikov kernel.
%
%   Example: Filter a microelectrode recording.
%
%      load MER.mat;
%      N       = length(x);
%      fss     = 750;                      % New sample rate (Hz)
%      Ns      = ceil(N*fss/fs);           % Number of samples
%      sti     = floor((fss/fs)*(si-1))+1; % Indices of spikes 
%      ks      = 1:Ns;                     % Sample index
%      xs      = zeros(Ns,1);              % Allocate memory for spike train
%      xs(sti) = 1;                        % Create spike train
%      KernelFilter(xs,fss,300);
%
%   J. McNames, "Optimal rate filters," EMBC 2005.
%
%   Version 1.00 JM
%
%   See also Lowpass, SmoothSeries, and RateFilter.

%====================================================================
% Process function arguments
%====================================================================
if nargin<1 | nargin>8,
    help KernelFilter;
    return;
    end;

fo = 200;                                                  % Default filter order (samples)
if exist('foa') & ~isempty(foa),
    fo = ceil(foa/2)*2;
    end;
  
pf = 0;                                % Default - no plotting
if nargout==0,                         % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Process Inputs
%====================================================================
x  = x(:);
nx = length(x);

%====================================================================
% Author-Specified Parameters
%====================================================================
sd = linspace(0.35,0.50,100)*(fs/fc);

%====================================================================
% Preprocessing
%====================================================================
ih = -fo/2:fo/2;                                           % Filter impulse response indices
nh = length(ih);
p2 = max(14,nextpow2(nh)+1);                               % Power of 2 to use for the FFT
nf = 2^(p2-1)+1;                                           % Number of useful frequencies from FFT
k  = (1:nf)';                                              % Indices of useful range of FFTs
f  = fs*(k-1)/2^p2;                                        % Frequencies
Hi = f<fc;                                                 % Ideal lowpass filter impulse response

hr = sinc(ih*2*fc/fs);                                     % The raw sinc function (natural scale)
hr = hr/sum(hr);

%====================================================================
% Optimize the Filter Width (Scale)
%====================================================================
ns = length(sd);
er = zeros(ns,1);
cf = 0;                                                    % Convergence flag

while ~cf,
    for c1=1:length(sd),
        %h  = exp(-ih.^2/(2*sd(c1).^2));                        % Gaussian impulse response
        h = (1-(ih/sd(c1)).^2).*(abs(ih)<sd(c1));
        h = h/sum(h);
        if h(1)~=0 || h(end)~=0,
            error('Filter order is insufficient.');
            end;
        H = fft(h(fo/2+1:end),2^p2);    
        H = real(H + conj(H))' - h(fo/2+1);
        H = H(k);
        er(c1) = mean(abs(Hi-H).^2);
        end;
    [ero,imin] = min(er);
    sdo = sd(imin);
    cf = 1;
    end;
    
%====================================================================
% Calculate the Filter's Impulse Response
%====================================================================    
h  = (1-(ih/sdo).^2).*(abs(ih)<sdo);
ho = h/sum(h);
        
%====================================================================
%Apply the Filter
%====================================================================
mx = mean(x);
y = filter(ho,1,[x-mx;zeros(nh,1)])+mx;                      
y = y(fo/2+(1:nx));

%====================================================================
% Calculate the Frequency Response and Stats
%====================================================================
if pf>=1 || nargout>2,
    fmn = 0;
    fmx = fs/2;

    p2 = max(14,nextpow2(nh)+1);                           % Power of 2 to use for the FFT
    nf = 2^(p2-1)+1;                                       % Number of useful frequencies from FFT
    k = (1:nf)';                                           % Indices of useful range of FFTs
    f = fs*(k-1)/2^p2;                                     % Frequencies

    H  = fft(ho(fo/2+1:end),2^p2);
    H  = H + conj(H) - h(1);
    H  = real(H);
    H  = H(k);
    H  = H(:);
    Hi = f<fc;                                             % The ideal frequency response 
    
    s = struct('DC',nan,'MSE',nan,'MSEPB',nan,'MSESB',nan,'MXEPB',nan,'MXESB',nan);
        
    pb0 = min(find(H<0.75)-1);                             % Magnitude-defined passband edge    
    sb0 = min(nf,max(find(H>0.25))+1);                     % Magnitude-defined stopband edge
    
    if isempty(pb0),
        pb0 = 1;
        end;
    
    ipb = 1:pb0;                                           % Passband indices
    isb = sb0:nf;                                          % Stopband indices
    itb = (ipb(end)+1):(isb(1)-1);                         % Transition band indices
    
    s.DC    = 2*sum(h)-h(1);                               % DC gain
    
    s.MSE   = mean((Hi-H).^2);                             % Total    mean squared error
    s.MSEPB = mean((1-H(ipb)).^2);                         % Passband mean squared error
    s.MSESB = mean((0-H(isb)).^2);                         % Stopband mean squared error
    
    s.MXEPB = max(abs(1-H(ipb)));                          % Passband maximum error
    s.MXESB = max(abs(0-H(isb)));                          % Stopband maximum error
    end;
    
%====================================================================
% Generate Plots of the Results
%====================================================================
if pf>=1,
    %----------------------------------------------------------------
    % Filtered Signal
    %----------------------------------------------------------------    
    figure;
    FigureSet(1);
    subplot(2,1,1);
        n = 1:nx;
        t = (n-0.5)/fs;
        ymx = max(y);
        ymn = min(y);
        yrg = ymx-ymn;
        ph  = stem(t,0.2*ymx*x);
        set(ph,'LineWidth',0.1);
        set(ph,'Marker','None');
        set(ph,'Color',[0.5 1.0 0.5]);    
        hold on;
            ph = plot(t,y,'k');
            set(ph,'LineWidth',1.5);
            hold off;
        xlim([0 nx/fs]);
        ylim([ymn-0.025*yrg ymx+0.025*ymx]);
        xlabel('Time (s)');
    subplot(2,1,2);
        h = plot(sd*fc/fs,er,'r');
        set(h,'LineWidth',2);
        hold on;
            h = plot(sdo*fc/fs,ero,'g.');
            set(h,'MarkerSize',25);
            hold off;
        xlim([sd(1) sd(end)]*fc/fs);
        xlabel('Filter Width (s)');
        ylabel('Mean Squared Error (MSE)');
    AxisSet;
        
    %----------------------------------------------------------------
    % Impulse & Frequency Response
    %----------------------------------------------------------------    
    figure;
    FigureSet(2);
    subplot(2,1,1);
        yrg = max(hr) - min(hr);
        ymn = min(hr) - 0.025*yrg;
        ymx = max(hr) + 0.025*yrg;
        tmx = nh/fs;
        
        t = ih/fs;
        hold on;        
            ph = plot(t,hr);
            set(ph,'LineWidth',3);
            set(ph,'Color',[0.5 1.0 0.5]);
            
            ph = plot(t,ho,'k');
            set(ph,'LineWidth',1.5);
            hold off;
        xlim([t(1) t(end)]);
        ylim([min(hr)-0.025*yrg max(hr)+0.025*yrg]);
        AxisLines;
        xlabel('Time (s)');
        ylabel('h(n)');    
        box off;
    subplot(2,1,2);       
        yrg = (1+s.MXEPB) - (0-s.MXESB);
        ymx = (1+s.MXEPB) + 0.05*yrg;
        ymn = (0-s.MXESB) - 0.05*yrg;
    
        ph = patch([0 f(ipb(end)) f(ipb(end)) 0],[1-s.MXEPB 1-s.MXEPB ymn ymn],'k');
        set(ph,'FaceColor',0.8*[1 1 1]);       
        set(ph,'LineStyle','None');       
    
       
        hold on;
            ph = patch([0 f(ipb(end)) f(ipb(end)) 0],[1+s.MXEPB 1+s.MXEPB ymx ymx],'k');
            set(ph,'FaceColor',0.8*[1 1 1]);       
            set(ph,'LineStyle','None');       
        
            ph = patch([f(isb(1)) f(end) f(end) f(isb(1))],[s.MXESB s.MXESB ymx ymx],'k');
            set(ph,'FaceColor',0.8*[1 1 1]);       
            set(ph,'LineStyle','None');          
           
            ph = patch([f(isb(1)) f(end) f(end) f(isb(1))],[-s.MXESB -s.MXESB ymn ymn],'k');
            set(ph,'FaceColor',0.8*[1 1 1]);       
            set(ph,'LineStyle','None');          
            
            ph = plot(f,f<fc,'k');
            set(ph,'LineWidth',1.0);
            
            ph = plot(f(ipb),H(ipb),'g',f(itb),H(itb),'y',f(isb),H(isb),'r');
            set(ph,'LineWidth',1.5); 
            
            %ph = plot([0 fs/2],b(1)*[1 1],'r',[0 fs],-b(2)*[1 1],'r');
            %set(ph,'LineWidth',0.5);
            
            hold off;
        xlim([0 fs/2]);
        ylim([ymn ymx]);
        AxisLines;
        xlabel('Frequency (Hz)');
        ylabel('|H(e^{j\omega})|');
        box off;
        set(gca,'Layer','Top');
    AxisSet;
    end;

%====================================================================
% Post Processing
%====================================================================
h = ho;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('y','h','s');
    end;