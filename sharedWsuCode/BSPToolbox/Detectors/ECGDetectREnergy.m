function [di] = ECGDetectREnergy(x,fsa,hra,pfa);
%DetectQRSEnergyBased: Detects QRS components of ECG signals
%
%   [di] = ECGDetectqREnergy(x,fs,hr,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=500 Hz.
%   hr   Expected range of the heart rate (Hz). Default=[0.8 3.4] Hz.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   di   Matrix of times (samples) of R components
%
%   Identifies the R time and amplitude for each individual beat. 
%   The amplitudes are measured relative to the detrended signal.
%   This is a complicated algorithm that attempts to locate the QRS
%   complexes based on amplitude and then corrects the times of the
%   estimated complexes using inter-beat variability.
%
%   Note that the times and amplitudes are returned because the
%   signal is internally resampled at approximately 1000 Hz. 
%
%   Example: Estimate ECG component times and amplitudes of an 
%   ECG signal.
%
%      load NoisyECG.mat; 
%      [di] = ECGDetectREnergy(ecg,fs,[],1);
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.408-412, 1996.
%
%   Version 1.00 JM
%
%   See also Detectors and HarmonicPSD.

%====================================================================
% Process Function Arguments
%====================================================================
if nargin<1,
    help ECGDetectREnergy;
    return;
    end;

nx = length(x);
if nx==0,
    fprintf('ERROR: Signal is empty.\n');
    return;
    end;

fs = 500;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;    
    
hr = [0.8 3.4];
if exist('hra') & ~isempty(hra) & length(hra)==2,
    hr = hra;
    end;    
    
pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Author-Specified Parameters
%====================================================================
Tvar   =  0.5;   % 0.3;  % Maximum expected beat-to-beat variation from median (percent)
EPD    =  0.70;  % Energy Peak Window length (seconds)

dwb    =   6;    % No. of seconds buffered before detection window
dwd    =  60;    % No. of seconds in detection window
dwa    =   6;    % No. of seconds buffered after detection window

wd     =   5;    % Window duration for spectral hear rate estimation (seconds)
nh     =  10;    % No. of harmonics to use in harmonic PSD estimation of heart rate

PRW    = 0.1;    % Maximum distance from energy peak to R Peak (seconds)
flp    = 2;      % Highpass drift Filter cutoff frequency (Hertz)
fhp    = 20;     % Lowpass noise Filter (hertz)

%====================================================================
% Preprocessing
%====================================================================
x     = x(:);                          % Make into a column vector
nx    = length(x);                     % No. Points in x
nm    = ceil(nx/(fs*60));              % Approximate No. minutes in x

nb    = ceil((nx/fs)*hr(2));           % Upper bound on number of beats
di    = zeros(nb,1);                   % Allocate memory for output matrix 

bimin = (fs)/hr(2);                    % Minimum expected samples per beat interval
bimax = (fs)/hr(1);                    % Maximum expected samples per beat

xmax  = max(x);                        % Maximum of input signal
xmin  = min(x);                        % Minimum of input signal

cb = 0;                                % Count of beats detected so far 
t = zeros(round(1.1*nm*60*max(hr)),1); % Allocate memory for output vector
r = zeros(round(1.1*nm*60*max(hr)),1); % Allocate memory for output vector

%====================================================================
% Main Loop - Calculate Raw Statistics
%====================================================================
for c1 = 1:nm, % Process minute by minute
	fprintf('Processing Minute %3d of %3d\n',c1,nm);
   
	%------------------------------------------------------------------
	% Check for clipping
	%------------------------------------------------------------------
	ki  = dwd*fs*(c1-1)+1;      % Windowed Segment initial index
	kf  = min([nx,dwd*c1*fs]);  % Windowed Segment final index    
	pr  = 0.1*(dwd*fs);         % Number of samples in 10% of the segment

    %------------------------------------------------------------------
	% Calculate statistics
	%------------------------------------------------------------------
	xmean = mean(x(ki:kf));       % Mean Value
	xstd  = std (x(ki:kf));       % Standard Deviation
   
	%------------------------------------------------------------------
	% Get indices and times of extended segment
	%------------------------------------------------------------------
	ki  = max([ 1,ki-round(dwb*fs)]); % Extended Segment initial index
	kf  = min([nx,kf+round(dwa*fs)]); % Extended Segment final index
   
	ti  = (ki-1)/fs;                  % Initial time of extended segment (0 based)
	tf  = (kf-1)/fs;                  % Final   time of extended segment 
   
	%dur = (kf-ki+1)/(dwd*fs);        % Duration (in minutes) of extended segment
    
    y = x(ki:kf);                     % Stage output
    
	%------------------------------------------------------------------
	% Remove Drift - Highpass Filter
	%------------------------------------------------------------------
	u = y - Lowpass(y,fs,flp,4);
    y = u;                            % Stage output
  
	%------------------------------------------------------------------
	% Impulse Filter - nonlinear squashing function
	%------------------------------------------------------------------
	%p   = percentile(y,[0.005 0.995]);
    p   = prctile(y,[0.5 99.5]);
	fc  = (p(2)+p(1))/2;             % Filter Center
	fr  = (p(2)-p(1));               % Filter Range
	if fr==0,
		fr = 1;
		end;
	u = fr*tanh((y-fc)/fr) + fc;      % Impulse Filter returns Signal, no impulse
    y = u;                            % Stage output
    
   	%------------------------------------------------------------------
	% Estimate Heart Rate 
	%------------------------------------------------------------------
    fmax = hr(2)*nh*2;                % Maximum signal frequency needed
    dsr  = floor(fs/fmax);            % Downsample rate
    u    = decimate(y,dsr);           % Decimate signal   
    wl   = round(wd*fs/dsr);          % Window duration

    %[p,f]    = pwelch(u,blackman(wl),[],2^15,fs/dsr);
    %[p,f]    = BlackmanTukey(u,fs/dsr,20,2^15);
    [p,f]    = Welch(u,fs/dsr,wd,[],2^15);
    %[hp,fh]  = HarmonicPSD(p,f,nh);
    hp = p;
    fh = f;
    
    id  = find(fh>=hr(1) & fh<=hr(2)); % Only consider HR in specified range
    fh  = fh(id);
    hp  = hp(id);
    
    % Find all peaks in the harmonic PSD
    k   = 2:length(hp)-1;
    fid = k((hp(k)> hp(k-1) & hp(k)>=hp(k+1)) | (hp(k)>=hp(k-1) & hp(k)> hp(k+1)));  
    fid = fid(find(hp(fid)>0.01*max(hp(fid))));  % Only consider those with power within 1% of max
    hrc = fh(fid); % Candidate heart rates
    
    % Find all peaks in the segment
    ny = length(y);
    k  = 2:ny-1;
    pi = k((y(k)>y(k-1) & y(k)>=y(k+1)) | (y(k)>=y(k-1) & y(k)>y(k+1)));  
    ps = -sort(-y(pi));    % Sort the peaks in ascending order
    pv = zeros(size(hrc)); % Peak variability
    for c2 = 1:length(hrc),
        enp    = round(hrc(c2)*(ny/fs)); % Expected number of peaks
        pv(c2) = var(ps(1:enp));
        end;
    [jnk,id] = min(pv);               % Find the candidate with the smallest peak amplitude variability
    ehr = hrc(id);                    % Estimated heart rate                  
    
    if 0,
        figure;
        FigureSet(1);
        subplot(2,1,1);
        h = plot(f,p);
        ylabel('PSD');
        subplot(2,1,2);
        h = plot(fh,hp,'b',fh(fid),hp(fid),'r.');
        ylabel('Harmonic PSD');
        xlabel('Frequency (Hz)');
        title(sprintf('Estimated HR: %5.3f Hz',ehr));                
        fprintf('Pausing...\n');
        pause;
        end;
        
    %------------------------------------------------------------------
	% Remove High-Frequency Noise - Lowpass Filter
	%------------------------------------------------------------------
	u = Lowpass(y,fs,fhp,4);
    y = u;                            % Stage output
    
	%------------------------------------------------------------------
	% Normalize Segment
	%------------------------------------------------------------------
	ymean = mean(y);
	ystd  = std(y);
	u     = (y - ymean)/ystd;         % Segment to be processed - normalized
    y     = u;                        % Stage output        
         
	%------------------------------------------------------------------
	% Find peaks in signal power 
	%------------------------------------------------------------------    
	fc = ehr;
    [pi,sp] = PowerPeaks(y,fs,fc,4,0);

    cnt = 1;
    ny  = length(y);                  % No. of samples in segment
    dur = ny/fs;                      % Duration of segment
    enp = ehr*dur;                    % Expected number of peaks
	np  = length(pi);                 % No. of peaks
    while ((np>1.5*enp) | (np<0.5*enp)) & cnt<3, % Do a dumb optimization on the cutoff frequency
		if np>1.5*enp,
            fc = fc/1.2;
        else
            fc = fc*1.2;
            end;
        [pi,sp] = PowerPeaks(y,fs,fc,4,0);
        np  = length(pi);    
		cnt = cnt + 1;
		end;

	if 0 & pf==2,
		fprintf('Number of Stage 1 loops: %d\n',cnt);        
		figure;
		FigureSet;
		k = 1:length(y);
        t = (k-1)/fs;
		plot(t,y,'g',t(pi),y(pi),'r.','MarkerSize',20);
		zoom on;
		fprintf('Pausing...\n'); 
        pause;
		end;       
        
	%------------------------------------------------------------------
	% Convert Segment Peak Indices to Signal Indices
	%------------------------------------------------------------------   
    t0 = dwd*(c1-1);       % Start boundary (sec) 
    t1 = dwd*c1;           % End boundary   (sec)
    k0 = t0*fs+1;          % Start boundary (samples)
    k1 = t1*fs+1;          % End boundary   (samples)
    ri = (pi-1) + ki;      % Translate segment indices into signal indices (samples)
    pi = ri(ri>k0 & ri<=k1); % Valid peak indices - those within the specified segment (samples)
    np = length(pi);       
    
	%------------------------------------------------------------------
	% Move Energy Peaks to Nearest Signal Peak
	%------------------------------------------------------------------      
    for c2 = 1:np,
        i0 = max( 1,pi(c2)-round(bimin/4)); % Index before peak
        i1 = min(nx,pi(c2)+round(bimin/4)); % Index after peak
        [mx,imx] = max(x(i0:i1));           % Find peak in signal
        if imx~=i0 & imx~=i1,               % Max should not be at an edge
            pi(c2) = (imx-1)+i0;
        else
            fprintf('WARNING: Peak %d found maximum at edge of boundary %d - %d\n',c2,i0,i1);
            end;
        end;
            
	%------------------------------------------------------------------
	% Assign to Final Set of Indices
	%------------------------------------------------------------------              
    di(cb+(1:np)) = pi;    % Times of detected peaks 
    cb = cb + np;          % Update count of No. of detected peaks    
    
    if 0,
        k = ki:kf;
        FigureSet(1,'wide');
        clf;
        id = di(1:cb);
        h = plot(k,x(k),'b',pi,x(pi),'r.');
        set(h(1),'Marker','.');
        set(h(1),'MarkerSize',15);
        set(h(2),'MarkerSize',6);
        xlim([ki kf]);
        zoom on;
        fprintf('Pausing...\n');
        pause;
        end;
	end;
    
%====================================================================
% Postprocessing
%====================================================================     
di = di(1:cb);

%====================================================================
% Plot Results
%====================================================================  
if pf~=0,    
    figure;
    FigureSet(1);
    subplot(2,1,1);
        tx = (0:nx-1)'/fs;
        plot(tx,x,'b',tx(di),x(di),'r.');
        xlim([0 (nx-1)/fs]);
        ylim([min(x) max(x)]);
        box off;
        ylabel('Signal (scaled)');
    subplot(2,1,2);
        np = length(di);
        ti = (di(1:np-1)+di(2:np))/2;
        ti = (ti-1)/fs;
        ib = diff(di)/fs;
        plot(ti,ib,'k.');
        box off;
        ylim([0 2]);
        xlim([0 (nx-1)/fs]);
        xlabel('Time (sec)');
        ylabel('Interbeat Intervals (s)');
    AxisSet;
    end;
        
%====================================================================
% Process Output Arguments
%====================================================================      
if nargout==0,
    clear('r','t');
    return;
    end;    