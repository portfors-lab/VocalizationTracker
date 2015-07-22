function [t,r] = ECGDetectREnergy(x,fsa,hra,pfa);
%DetectQRSEnergyBased: Detects QRS components of ECG signals
%
%   [t,a] = ECGDetectqREnergy(x,fs,hr,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=500 Hz.
%   hr   Expected range of the heart rate (Hz). Default=[0.8 3.4] Hz.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   t    Matrix of times (sec) of R components
%   r    Matrix of amplitudes  of R components 
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
%      load NOISYECG.mat; 
%      [t,a] = ECGDetectREnergy(noisyecg,fs,[],1);
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.408-412, 1996.
%
%   Version 1.00 JM
%
%   See also Detectors and HarmonicPSD.

if IsScript(mfilename),
    clear all;
    load NOISYECG;
    x = noisyecg(550e3:600e3);
    fsa = fs;
    pfa = 2;
    end;

%====================================================================
% Process Function Arguments
%====================================================================
if ~IsScript(mfilename) & nargin<1,
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
if IsScript(mfilename) | nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Author-Specified Parameters
%====================================================================
fst    = 1000;   % Target sample rate
Tvar   =  0.5;   % 0.3;  % Maximum expected beat-to-beat variation from median (percent)
EPD    =  0.70;  % Energy Peak Window length (seconds)

dwb    =   6;    % No. of seconds buffered before detection window
dwd    =  60;    % No. of seconds in detection window
dwa    =   6;    % No. of seconds buffered after detection window
wd     =   8;    % Window duration for spectral hear rate estimation (seconds)

PRW    = 0.1;    % Maximum distance from energy peak to R Peak (seconds)
flp    = 2;      % Highpass drift Filter cutoff frequency (hertz)
fhp    = 20;     % Lowpass noise Filter (hertz)

%====================================================================
% Calculate Initial Values & Constants
%====================================================================
usr    = ceil(fst/fs);               % Upsample rate
fse    = usr*fs;                     % Effective sample rate
nx     = length(x);                  % No. Points in x
nm     = ceil(nx/(fs*60));           % Approximate No. minutes in x

nb     = ceil((nx/fs)*(1/60)*hr(2)); % Upper bound on number of beats
T      = zeros(nb,5);                % Allocate memory for output matrix 
A      = zeros(nb,5);                % Allocate memory for output matrix

bimin = (usr*fs)/hr(2);              % Minimum expected samples per beat interval
bimax = (usr*fs)/hr(1);              % Maximum expected samples per beat

xmax   = max(x);                     % Maximum of input signal
xmin   = min(x);                     % Minimum of input signal

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
	% Upsample Segment
	%------------------------------------------------------------------
    u   = interp(y,usr);              % Upsample the segment
	u   = u(1:length(u)-(usr-1));     % Drop the last (usr-1) samples
	kiu = (ki-1)*usr+1;               % Recalculate effective initial index
	kfu = (kf-1)*usr+1;               % Recalculate effective final index 
    yu  = u;                          % Save stage output in yu for last stage
    y   = u;                          % Stage output
  
	%------------------------------------------------------------------
	% Impulse Filter - nonlinear squashing function
	%------------------------------------------------------------------
	p    = percentile(y,[0.005 0.995]);
	fc   = (p(2)+p(1))/2;             % Filter Center
	fr   = (p(2)-p(1));               % Filter Range
	if fr==0,
		fr = 1;
		end;
	u = fr*tanh((y-fc)/fr) + fc;      % Impulse Filter returns Signal, no impulse
    y = u;                            % Stage output
    
   	%------------------------------------------------------------------
	% Estimate Heart Rate 
	%------------------------------------------------------------------
    fmax = hr(2)*10*2;                % Maximum signal frequency needed
    dsr  = floor(fse/fmax);           % Downsample rate
    u    = decimate(y,dsr);           % Decimate signal   
    wl   = round(wd*fse/dsr);         % Window duration

    [p,f]    = pwelch(u,blackman(wl),[],2^15,fse/dsr);
    [hp,fh]  = HarmonicPSD(p,f,10);
    
    id  = find(fh>=hr(1) & fh<=hr(2)); % Only consider HR in specified range
    fh  = fh(id);
    hp  = hp(id);
    [jnk,id] = max(hp);               % Find maximum    
    ehr = fh(id);                      % Estimated heart rate
    
    if id==1 | id==length(id),
        fprintf('WARNING: Estimated HR outside of expected range.\n');
        end;

    %------------------------------------------------------------------
	% Remove High-Frequency Noise - Lowpass Filter
	%------------------------------------------------------------------
	u = Lowpass(y,fse,fhp,4);
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
    [pi,sp] = PowerPeaks(y,fse,fc,4,0);

    cnt = 1;
    ny  = length(y);                  % No. of samples in segment
    dur = ny/fse;                     % Duration of segment
    enp = ehr*dur;                    % Expected number of peaks
	np  = length(pi);                 % No. of peaks
    while ((np>1.5*enp) | (np<0.5*enp)) & cnt<3, % Do a dumb optimization on the cutoff frequency
		if np>1.5*enp,
            fc = fc/1.2;
        else
            fc = fc*1.2;
            end;
        [pi,sp] = PowerPeaks(y,fse,fc,4,0);
        np  = length(pi);    
		cnt = cnt + 1;
		end;

	if 0 & pf==2,
		fprintf('Number of Stage 1 loops: %d\n',cnt);        
		figure;
		FigureSet;
		k = 1:length(y);
        t = (k-1)/fse;
		plot(t,y,'g',t(pi),y(pi),'r.','MarkerSize',20);
		zoom on;
		fprintf('Pausing...\n'); 
        pause;
		end;       
        
	%------------------------------------------------------------------
	% Move Energy Peaks to Nearest Signal Peak
	%------------------------------------------------------------------   
    ny = length(y);
    for c2 = 1:np,
        i0 = max( 1,pi(c2)-round(bimin/4)); % Index before peak
        i1 = min(ny,pi(c2)+round(bimin/4)); % Index after peak
        [mx,imx] = max(y(i0:i1));           % Find peak in signal
        if imx~=i0 & imx~=i1,               % Max should not be at an edge
            pi(c2) = (imx-1)+i0;
            end;
        end;

	%------------------------------------------------------------------
	% Update output vector
	%------------------------------------------------------------------   
    t0 = dwd*(c1-1);      % Start boundary 
    t1 = dwd*c1;          % End boundary
    tp = ((kiu:kfu)-1)/fse;      
    ti = tp(pi);
    id = find(ti>t0 & ti<=t1); % Indices of peaks within the interval under consideration
    pi = pi(id);
    np = length(id);
    
    t(cb+(1:np)) = ti(id);     % Times of detected peaks 
    r(cb+(1:np)) = yu(pi); % Amplitudes of detected peaks (after drift removal)
    cb = cb + np;          % Update count of No. of detected peaks    
	end;
    
%====================================================================
% Postprocessing
%====================================================================     
t = t(1:cb);
r = r(1:cb);

%====================================================================
% Plot Results
%====================================================================  
if pf~=0,    
    figure;
    FigureSet(1);
    subplot(2,1,1);
        tx = (0:nx-1)'/fs;
        plot(tx,x-mean(x),'b',t,r,'r.');
        xlim([0 (nx-1)/fs]);
        ylim([min(x-mean(x)) max(r)]);
        box off;
        ylabel('Signal (scaled)');
    subplot(2,1,2);
        np = length(t);
        ti = (t(1:np-1)+t(2:np))/2;
        ib = diff(t);
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
if ~IsScript(mfilename) & nargout==0,
    clear('r','t');
    return;
    end;    