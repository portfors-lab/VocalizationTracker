function [B,t,d,pa] = Beatogram(x,fsa,wla,doa,rga,nsa,pfa);
%Beatogram: Plots pulse morphology versus time as an image
%
%   [B,t,d,pa] = Beatogram(x,fs,do,rg,ns,pf);
%
%   x    Input signal       
%   fs   Sample rate in (Hz). Default = 125 Hz.
%   wl   Length of window to use (s). Default = 8 s. 
%   do   Order of derivative to show (integer 0-2). Default = 0.
%   rg   Range of each beat to show (s). Default = [-1 1].
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        Default = min(2^10,NX).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   B    Matrix containing the image created
%   t    Times (columns) at which B was estimated (s)
%   d    Delays (rows) at which B was estimated (s)
%   pa   Average number of peaks detected in each window
%
%   Generates an image that shows how the pulse morphology changes over
%   time. If the order is specified, then Beatogram plots the 1st or
%   2nd derivative of each pulse.
%
%   Example: Plot the Beatogram of an intra-cranial pressure (ICP)
%   signal.
%
%      load ICP;
%      Beatogram(icp,fs);
%
%   J. McNames, J. Bassale, M. Aboy, C. Crespo, B. Goldstein, "Techniques 
%   for the Visualization of Nonstationary Biomedical Signals," accepted 
%   for presentation at Biosignal 2002. 
%
%   Version 1.00 JM
%
%   See also Spectrogram and Cohereogram.

%====================================================================
% Process function arguments, fill in defaults
%====================================================================
if nargin<1,
    help Beatogram;
    return;
    end;

nx = length(x); % Length of input vector
if nx==0,
    fprintf('Error: Signal is empty.\n');
    return;
    end;

fs = 125; % Default sampling rate (Hz)
if exist('fsa') & ~isempty(fsa),
    fs = max(fsa,0);
    end;
    
wl = round(8*fs); % Default: 8 seconds in terms of samples
if exist('wla') & ~isempty(wla),
    wl = round(wla*fs);
    wl = max(150,wl); % Can't be less than 150 samples
    wl = min(wl,nx);  % Can't be larger than nx
    end;    

do = 0; % Derivative order. Default: none (0)
if exist('doa') & ~isempty(doa),
    do = round(doa);
    do = max(do,0);
    do = min(do,2);
    end;
    
rg = [-1 1]; % Default range: -1 to 1 sec
if exist('rga') & ~isempty(rga) & length(rga)==2,
    rg = rga;
    rg(1) = min(rg(1),0); % Should be non-positive
    rg(2) = max(rg(2),0); % Should be non-negative
    end;          
    
ns = min(500,nx); % Default: 500 pixels
if exist('nsa') & ~isempty(nsa),
    ns = max(10,min(nsa,nx));
    end;       
          
pf = 0;                                 % Default - no plotting
if nargout==0,                          % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing
%====================================================================    
if size(x,1)~=1, % Make into a row vector for the mean operator
    x = x';
    end;

xlp = Lowpass(x,fs,0.01*fs/2)'; % Remove trend
xdt = x-xlp;                  % Detrended

bi = PressureDetect(x,fs);    % Returns indices of detected beats
st = (length(x)-1)/(ns-1);    % Step size
st = max(st,1);               % Don't let step size be less than a sample
k  = 1:st:nx;                 % Evaluation times (samples)
ne = length(k);               % No. of evaluation times
ll = min(-1,floor(rg(1)*fs)); % Lower limit, no greater than -1 sample
ul = max( 1,ceil(rg(1)*fs));  % Upper limit, no less than 1 sample
pw = ll:ul;                   % Peak indices window (without offset)
nv = length(pw)-do;           % Number of vertical pixels
B  = zeros(nv,ns);            % Beatogram image matrix allocation

%====================================================================
% Main Loop
%====================================================================  
pa = 0; % Counter of number of peaks detected
for c = 1:ne,
    ke  = k(c);                                   % Evaluation time
    i0  = max( 1-min(pw),round(ke-wl/2));         % Lower window limit
    i1  = min(nx-max(pw),round(ke+wl/2));         % Upper window limit
    pi  = bi(bi>=i0 & bi<=i1);                    % Select peak indices
    np  = length(pi);                             % No. peaks found
    if np>0,
        ps  = xdt(pi*ones(1,nv) + ones(np,1)*pw); % Selected pulses     
        pm  = median(ps,1)';                      % Median pulse
    else                                          % If no peaks detected
        pm = zeros(nv,1);
        end;
    if do==1,
        pm = fs*diff(pm);                        % Velocity
    elseif do==2,
        pm = fs^2*diff(diff(pm));                % Acceleration
        end;
    B(:,c) = pm;
    pa     = pa + np;
    end;
pa = pa/ne;
%fprintf('Average of %5.1f peaks per slice.\n',npa);

%====================================================================
% Postprocessing
%====================================================================  
t = (k-1)/fs;

if do==0,
    d = pw/fs;
elseif do==1,
    np = length(pw);
    d  = (pw(1:np-1) + pw(2:np))/(2*fs);
elseif do==2,
    np = length(pw);
    d  = pw(2:np-1)/fs;
    end;
    
t = t(:); % Convert to column vector
d = d(:); % Convert to column vector

%====================================================================
% Plot Results
%====================================================================  
if pf==1 | IsScript(mfilename),    
    LB  = 0.08;   % Left boundary
    RB  = 0.98;   % Right boundary
    WD  = RB-LB;  % Width
    NSA = 5;      % No. segment axes
    NSP = 2*fs+1; % No. segment points in each axis - must be odd
    
    figure;
    FigureSet;
    
    % Top axis
    ha1 = axes('Position',[LB 0.33 WD 0.63]); % Beatogram
        p  = prctile(xdt,[2 25 75 99]);
        zmin = p(1) - 0.30*(p(3)-p(2));
        zmax = p(4) + 0.30*(p(3)-p(2));
        imagesc(t,d,B,[zmin zmax]);
        ylabel('Morphology (s)');
        title('Beatogram: Pulse Morphology versus Time');
        set(gca,'XTickLabel','');
        set(gca,'Box','Off');
        
    % Middle axis
    ha2 = axes('Position',[LB 0.22 WD 0.10]);
        kr = 1:nx; 
        tr = (kr-1)/fs; % Times of raw signal
        h = plot(tr, x);
        set(h,'Color',[0.7 0.7 1.0]); % Light blue
        hold on;
            h = plot(tr,xlp);              % Trend
            h = plot(tr(bi), x(bi), 'r.'); % Mark detected beats
            set(h,'MarkerSize',2);         % Small marks
            hold off;
        p = prctile(x,[2 25 75 99]);
        ymin2 = p(1) - 0.30*(p(3)-p(2));
        ymax2 = p(4) + 0.30*(p(3)-p(2));
        xlim([0 max(tr)]);
        ylim([ymin2 ymax2]);

    % Bottom axis: Plot derivatives in the subsegments
    if do==0,
        sg = xdt;
    elseif do==1,
        sg = diff(xdt)*fs;
    elseif do==2,
        sg = diff(diff(xdt))*fs^2;
        end;
    
    p = prctile(sg,[2 25 75 99]);
    ymin = p(1) - 0.30*(p(3)-p(2));
    ymax = p(4) + 0.30*(p(3)-p(2));
    for cnt = 1:NSA,
        bf = 0.01;                         % Buffer between plots 
        wd = ((RB-LB) - 0.01*(NSA-1))/NSA; % Axis width
        le = LB + (wd+bf)*(cnt-1);
        ha = axes('Position',[le 0.08 wd 0.10]);   
        cp = round(cnt*nx/(NSA+1));        % Center point of axis
        [jnk,I] = min(abs(cp-bi));         % Find the closest peak
        cp = bi(I);
        k  = cp + (-(NSP-1)/2:(NSP-1)/2);
        ts = (k-1)/fs;
        h  = plot([min(ts) max(ts)],[0 0],'k:');
        hold on;
            h = plot((cp-1)/fs*[1 1],[ymin ymax],'k:');
            if do==0,
                sg = xdt(k);
                h = plot(ts, sg,'b');
                h = plot(tr(bi), xdt(bi), 'r.');
            elseif do==1,
                np = length(ts);
                ts = (ts(1:np-1) + ts(2:np))/2;
                sg = diff(xdt(k))*fs;
                h = plot(ts, sg,'b');
            elseif do==2,
                np = length(ts);
                ts = ts(2:np-1);
                sg = diff(diff(xdt(k)))*fs^2;
                h = plot(ts, sg,'b');
                end;
            set(h,'MarkerSize',8);
            hold off;
        xlim([min(ts) max(ts)]);
        ylim([ymin ymax]);
        if cnt>1,
            set(gca,'YTickLabel','');
            end;
        if cnt==(NSA-1)/2+1,
            xlabel('Time (s)');
            end;
        axes(ha1);
        hold on;
            plot((cp-1)/fs*[1 1],[min(d) max(d)],'k:');
            hold off;
        axes(ha2);
        hold on;
            plot((cp-1)/fs*[1 1],[ymin2 ymax2],'k:');
            hold off;
        end;    
    axes(ha1);
    zoom on;
    AxisSet;
    end;

if ~IsScript(mfilename) & nargout==0,
    clear('B','t','d','pa');    
    end;
