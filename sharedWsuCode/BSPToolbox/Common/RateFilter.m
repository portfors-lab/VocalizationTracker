function [y,h,s] = RateFilter(x,fs,fc,foa,rpa,lba,mda,pfa)
%RateFilter: Rate estimating filter for point processes
%
%   [y,h,s] = RateFilter(x,fs,fc,fo,rp,lb,md,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz).   
%   fc   Cutoff frequency (Hz).
%   fo   Filter order (samples). Must be even. Default = 200.
%   rp   Maximum ripple (% of ideal response max ripple). 
%        Default = 25%.
%   lb   Impulse response lower bound (% of ideal response min). 
%        Default = 25%.
%   md   Maximum deviation of response at DC (mean) from 1 (%). 
%        Default = 1%.
%   pf   Plot flag: 0=none (default), 1=screen
%
%   y    Filtered signal (estimated instantaneous event rate).
%   h    Filter impulse response.
%   s    Structure containing filter statistics.
%
%   Creates an optimal FIR intensity filter that is as close as 
%   possible to an ideal lowpass filter subject to the ripple and 
%   negativity constraints.
%
%   Example: Filter a microelectrode recording.
%
%      load MER.mat;
%      N       = length(x);
%      fss     = 750;                      % New sample rate (Hz)
%      Ns      = ceil(N*fss/fs);           % Number of samples
%      sti     = floor((fss/fs)*(si-1))+1; % Indices of spikes 
%      ks      = 1:Ns;                     % Sample index
%      ts      = (ks-0.5)/fss;             % Time index
%      xs      = zeros(Ns,1);              % Allocate memory for spike train
%      xs(sti) = 1;                        % Create spike train
%      RateFilter(xs,fss,50);
%
%   J. McNames, "Optimal rate filters," EMBC 2005, in review.
%
%   Version 1.00 JM
%
%   See also Lowpass, SmoothSeries, and KernelFilter.

%====================================================================
% Process function arguments
%====================================================================
if nargin<1 | nargin>8,
    help RateFilter;
    return;
    end;

fo = 200;                                                  % Default filter order (samples)
if exist('foa') & ~isempty(foa),
    fo = ceil(foa/2)*2;
    end;
 
rp = 0.25;                                                 % Maximum ripple
if exist('rpa') & ~isempty(rpa),
    rp = rpa/100;
    end;
    
lb = 0.25;                                                 % Minimum possible impulse response (% of maximum)
if exist('lba') & ~isempty(lba),
    lb = lba/100;
    end;
    
md = 0.01;                                                 % Maximum possible deviation of the average
if exist('mda') & ~isempty(mda),
    md = mda/100;
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
% Preprocessing
%====================================================================
nh = fo/2+1;                                               % Length of one-sided impulse response 
n  = (0:nh-1)';                                            % Sample indices of the impulse response
hr = sinc(n*2*fc/fs);                                      % The raw sinc function (natural scale)
hs = 2*sum(hr)-hr(1);                                      % Two-sided sum of the impulse response
hr = hr/hs;                                                % Raw, scaled sinc function

hr1    = hr;                                               % Create version with first element scaled 
hr1(1) = sqrt(1/2)*hr(1);                                  % to account for the fact that it only shows up once

%====================================================================
% Optimization Parameters
%====================================================================
op     = optimset('LargeScale','Off','Display','off');
C      = eye(nh);
C(1,1) = sqrt(1/2);
nc     = nh*(ceil(2/fc));
A      = zeros(nc,nh);
b      = zeros(nc,1);

cr      = 1;
A(cr,:) = 2*ones(1,nh);
A(cr,1) = 1;
b(cr,:) = 1+md;

cr      = 2;
A(cr,:) = -2*ones(1,nh);
A(cr,1) = -1;
b(cr,:) = -(1-md);

nc   = cr;
A    = A(1:nc,:);   
b    = b(1:nc);

hlb     = -inf*ones(nh,1);                                 % Allocate memory & initialize lower bound on h
hub     =  inf*ones(nh,1);                                 % Allocate memory & initialize upper bound on h
ib      = find(n>=0.50*fs/fc);                              % Values of n where the bounds apply (outside of the main lobe)
iub     = find(n>=0.75*fs/fc);
irb     = find(n>=1.00*fs/fc);
hrmin   = min(hr);                                         % Minimum value of ideal impulse response (also the maximum ripple)
if min(lb,rp)>=1.0,
    hlb = [];
else
    hlb      = lb*hrmin*ones(nh,1);                         % Lower bound of h
    hlb(irb) = min(lb,rp)*hrmin;                            % Lower bound of h with ripple
    end;
if rp>=1.0,
    hub = [];
else
    hub(iub) = rp*abs(hrmin);                                   % Upper bound on h
    end;
    
%====================================================================
% Solve for Optimal Filter
%====================================================================
hi     = hr;                                               % Use the optimal filter response as a starting point 
if ~isempty(hlb),
    id     = find(hi<hlb);                                     % Find the points where hi < lower bound
    hi(id) = hlb(id);                                          % Clip the initial impulse response from below
    end;
if ~isempty(hub),
    id     = find(hi>hub);                                     % Find the points where hi > upper bound
    hi(id) = hub(id);                                          % Clip the initial impulse response from above   
    end;

[ho,rn,rs,ef,om] = lsqlin(C,hr1,A,b,[],[],hlb,hub,hi,op);

if ef~=1,
    switch ef,
    case 3,
        warning('Change in the residual smaller that the specified tolerance.');
    case 0,
        warning('Maximum number of iterations exceeded.');
    case -2,
        warning('Problem is infeasible.');
    case -4,
        warning('Ill-conditioning prevents further optimization.');
    case -7,
        warning('Magnitude of search direction became too small; no further progress can be made.');
    otherwise
        warning(om.message);
        end;
    end;
        
%====================================================================
%Apply the Filter
%====================================================================
h = [flipud(ho(2:end));ho];
y = filter(h,1,[x;zeros(nh,1)]);                                         % Still need to deal with end conditions
y = y((nh-1)+(1:nx));
y = y*fs;

%====================================================================
% Calculate the Frequency Response and States
%====================================================================
if pf>=1 || nargout>2,
    fmn = 0;
    fmx = fs/2;

    p2 = max(14,nextpow2(nh)+1);                           % Power of 2 to use for the FFT
    nf = 2^(p2-1)+1;                                       % Number of useful frequencies from FFT
    k = (1:nf)';                                           % Indices of useful range of FFTs
    f = fs*(k-1)/2^p2;                                     % Frequencies

    H  = fft(ho,2^p2);
    H  = H + conj(H) - ho(1);
    H  = real(H);
    H  = H(k);
    Hi = f<fc;                                             % The ideal frequency response 
    
    s = struct('DC',nan,'MSE',nan,'MSEPB',nan,'MSESB',nan,'MXEPB',nan,'MXESB',nan);
    
    
    pb0 = min(find(H<0.75)-1);                             % Magnitude-defined passband edge
    %pb0 = min(find(H<2-max(H(1:pb0)))-1);    
    
    sb0 = max(find(H>0.25))+1;                             % Magnitude-defined stopband edge
    %sb0 = max(find(H>abs(min(H(sb0:end)))))+1; 
    
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
        
        n = 0:nh-1;
        t = n(:)/fs;
        ph = patch([min(t(ib)) tmx tmx min(t(ib))],[max(hlb) max(hlb) ymn ymn],'k');
        set(ph,'FaceColor',0.8*[1 1 1]);       
        set(ph,'LineStyle','None');       
        hold on;
            ph = patch([min(t(ib)) tmx tmx min(t(ib))],[min(hub) min(hub) ymx ymx],'k');
            set(ph,'FaceColor',0.8*[1 1 1]);       
            set(ph,'LineStyle','None');        
        
            ph = plot(t,hr);
            set(ph,'LineWidth',3);
            set(ph,'Color',[0.5 1.0 0.5]);
            
            ph = plot(t,ho,'k');
            set(ph,'LineWidth',1.5);
            hold off;
        xlim([0 tmx]);
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
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('y','h','s');
    end;
    
% function [f,g,H] = myfun(h,hr)
%     f = sum((h-hr).^2);
%     g = 2*(h-hr);    
%     H = 2*diag(h);
%     disp(nargout);
