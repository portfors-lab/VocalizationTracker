function [ST] = APDetect(x,fs,mr,pf);
%APDetect: Detects extracellular action potentials.
%
%   [a] = APDetect(x,fs,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=20,00 Hz.
%   mr   Maximum expected firing rate (Hz). Default = 100.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   a    Vector of detected spike times (samples)
%
%   Identifies the times at which action potentials occur. The times
%   identified are the times at which the action potentials reach
%   their maximum.
%
%   Example: Estimate ECG component times and amplitudes of an 
%   ECG signal.
%
%      load Tremor.mat; 
%      [a] = APDetect(x,fs);
%
%   J. McNames, "Microelectrode Signal Analysis Techniques for 
%   Improved Localization," submitted to Microelectrode Recordings
%   in Movement Disorder Surgery, Theime, 2002.
%
%   Version 2.00.22 JM
%
%   See also PowerPeaks and Detectors.

%close all;
%clear all;
%error('This function is under development. JM Jul 02');

%====================================================================
% Process Function Arguments
%====================================================================
if nargin<1,
    help APDetect;
    return;
    end;

nx = length(x);
if nx==0,
    error('Signal is empty.');
    end;

fs = 20000;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end; 
    
mr = 100;
if exist('mra') & ~isempty(mra),
    mr = mra;
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

%====================================================================
% Preprocessing
%====================================================================
x = x(:);                          % Convert to column vector
x = (x-mean(x))/std(x);            % Normalize to have zero mean and unit standard deviation

%====================================================================
% Main Loop
%====================================================================
[pi,sp] = PowerPeaks(x,fs,mr);

[B,BD,Step] = FindModes(sp,ATH,NPP,MTH);
nb  = size(B,1);                   % Number of Bumps
if nb>=2,
	ETh = B(nb-1,3);               % Set Energy Threshold
else
	ETh = B(1,3);
	end;
ESI = PI(EP>ETh);                  % Energy Spike Indices

% Eliminate spikes that are too close together (<TSM sec apart)
sm  = (TSM/1000)*FS;
dsi = diff(ESI);
k   = 1:length(dsi);
k   = k(dsi<sm);
for c1 = 1:length(k),
	i1 = ESI(k(c1)  );
	i2 = ESI(k(c1)+1);
	if SE(i1)>SE(i2),
		ESI(k(c1)+1) = -1;
	else
		ESI(k(c1)) = -1;
		end;
	end;
if ~isempty(ESI),
	ESI  = ESI(ESI~=-1);	
	end;

NES = length(ESI);                 % Number of Energy Spikes

fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);

if NES<NSM,
	fprintf('Breaking because of insufficient spikes detected.	\n');
	return;
	end;



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


PLOT = 1;
if exist('PlotArg') & ~isempty(PlotArg),
    PLOT = PlotArg;
    end;

PRINT = 1;      % Print plots?
if exist('PrintArg') & ~isempty(PrintArg),
   PRINT = PrintArg;
   end;
      
GREEN = [0.00 0.70 0.00];   % Color Vector (RGB) for green
BROWN = [0.75 0.47 0.09];   % Color Vector (RGB) for brown
BROWN = [0.38 0.24 0.05];   % Color Vector (RGB) for brown
     
FHP = 200;      % (Hz) Cutoff Frequency for High Pass Filter          

ATH = 0.005;    % 0.005 Minimum bump area (as proportion of total) for FindBumps
NPP = 4;        % 4     Maximum number of expected bumps
MTH = 0.75;     % 0.85  Minimum Threshold

EDR = 2.0;      % 2.0   (ms) Duration of window for energy location of spikes

NTL = 3;        % 3     Number of Template Loops in Phase II of Template Construction
OTP = 0.25;     % 0.25  Template offset (+ = move window to right of peak)
TDL = 1.0;      % 1.0   (ms) Approximate Duration of a Template, default
TSM = 1.0;      % (ms)  Minimum allowed separation between spikes 
NTS = 250;      % 250   Maximum number of spikes used to Construct Template

IOT = 0.001;    % (sec) Threshold to declare overlap 

%NTL = 10;       
%OTP = 0.75;     % Template Offset (+ = move window to right of peak)
%TDL = 0.0040;   %

NTP = 50;       % Maximum Number of points in a template

CFM = 0.1;      % Cutoff Frequency Multiplier for Estimating Spike Density
FTR = [4 7 10]; % (Hz) Range of Tremor Frequencies

NSM = 25;       % Minimum number of spikes required to proceed with analysis

PRT = 0.975;    % Pause Rate Threshold

%=========================================================
% Upsample the Signal
%=========================================================
NR = length(YR);

if FSR<25e3,
	USR = ceil(25e3/FSR);
	fprintf('Upsampling by %d...',USR);
	t0  = clock;
	Y   = interp(YR,USR);
	Y   = Y(1:USR*NR-1);
	FS  = FSR*USR;
	fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);
else
    FS  = FSR;
	Y   = YR;
	end;
Y  = Y - mean(Y);
NY = length(Y);

%=========================================================
% Raw & Upsampled Signal
%=========================================================
if 0 & PLOT==1,
    figure;
	FigureSet(1,FGWD,FGHT);
    k  = round(5*FS /FHP):round(20*FS /FHP); % 50 Time Constants
  	kr = round(5*FSR/FHP):round(20*FSR/FHP); % 50 Time Constants
	t  = (0:NY-1)/FS;
    tr = (0:NR-1)/FSR;
    h  = plot(tr(kr),YR(kr),'r',t(k),Y(k),'b');
    set(h,'LineWidth',1.2*LWM);
    set(h,'Marker','.');
    set(h(1),'MarkerSize',20);
	tmn = min(t(k));
	tmx = max(t(k));
	ymn = 1.1*min(min(Y(k)),min(YR(kr)));
	ymx = 1.1*max(max(Y(k)),max(YR(kr)));
	axis([tmn tmx ymn ymx]);
	xlabel('Time (sec)');
	ylabel('Raw Signal (Not Normalized)');
	st = sprintf(': Upsampling Sample %3d Hz   Before (red) After (blue)',FHP);
	title([NAME st]);
	set(gca,'Box','off');
	axisset(AFS);  
	zoom on;
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'Upsample');
		eval(st);
	elseif PRINT==2,
		print; 	
		end;
   end;
   
%=========================================================
% High-Pass Filter the Signal
%=========================================================
fprintf('High Pass Filtering...');
t0  = clock;
XHP = HighPass(  Y,FHP/(FS/2));
XHP = HighPass(XHP,FHP/(FS/2));
fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);

%=========================================================
% Plot Filtered & Unfiltered Signal
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);
	k = 1:round(200*FS/FHP); % 50 Time Constants
	k = 1:length(Y); % 50 Time Constants
	t = (0:NY-1)/FS;
	h = plot(t(k),Y(k),'r',t(k),XHP(k),'b');
	set(h,'LineWidth',1.2*LWM);
	tmn = min(t(k));
	tmx = max(t(k));
	ymn = 1.1*min(min(Y(k)),min(XHP(k)));
	ymx = 1.1*max(max(Y(k)),max(XHP(k)));
	axis([tmn tmx ymn ymx]);
	xlabel('Time (sec)');
	ylabel('Raw Signal (Not Normalized)');
	st = sprintf(': High Pass Filter Sample %3d Hz   Before (red) After (blue)',FHP);
	title([NAME st]);
	set(gca,'Box','off');
	axisset(AFS);  
	zoom on;
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'HighPass');
		eval(st);
	elseif PRINT==2,
		print; 	
		end;
   end;
     
%=========================================================
% Plot the Estimated Signal Spectrum
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	fprintf('Estimating power spectral density...');
	t0 = clock;

	ns = 2^floor(log2(NY));
	k  = 1:(ns/(2*128));
	f  = (FS/2)*(k-1)/length(k);

	Py = bart2(Y(1:ns),128,FS);
	Py = Py(k);
	%Py = Py*var(Y(1:ns))/(sum(Py)*(ns/FS));
	h  = plot(f,Py(k),'r');
	set(h,'LineWidth',1.8*LWM);

	hold on;
		%Px = bart2(XHP(1:ns),128,FS);
		Px = bart2(XHP(1:ns),128,FS);
		Px = Px(k);
      %Px = Px*var(XHP(1:ns))/(sum(Px)*(ns/FS));
      %Px(1:2) = 0;
		h  = plot(f,Px(k),'b');
		set(h,'LineWidth',1.5*LWM);
		hold off;

	fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);

	axis([0 FSR/2 0 1.05*max(Px(k))]);
	xlabel('Frequency (Hz)');
	ylabel('Power Spectral Density');
	st = sprintf(': Normalized PSD   Raw (Red) Filtered (Blue)   Cutoff:%3d Hz',FHP);
	title([NAME st]);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'SignalPSD');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
   end;
   
%=========================================================
% Normalize the signal
%=========================================================
X  = (XHP-mean(XHP))/std(XHP); 
[nr,nc] = size(X);
if nc==1,
   X = X';
	end;

NX = length(X);

%=========================================================
% Spike Energy
% Calculate Spike Energy, Find Peaks, and Pick Threshold
%=========================================================
fprintf('Finding energy spikes...');
t0 = clock;

EDR = EDR/1000;                    % Convert energy duration of energy window to ms
NEP = ceil(EDR*FS);                % Number of energy points
[PI,SE] = EnergySpikes(X,NEP);
EP = SE(PI);                       % Energy Peaks

[B,BD,Step] = FindBumps(EP,ATH,NPP,MTH);
nb  = size(B,1);                   % Number of Bumps
if nb>=2,
	ETh = B(nb-1,3);                 % Set Energy Threshold
else
	ETh = B(1,3);
	end;
ESI = PI(EP>ETh);                  % Energy Spike Indices

% Eliminate spikes that are too close together (<TSM sec apart)
sm  = (TSM/1000)*FS;
dsi = diff(ESI);
k   = 1:length(dsi);
k   = k(dsi<sm);
for c1 = 1:length(k),
	i1 = ESI(k(c1)  );
	i2 = ESI(k(c1)+1);
	if SE(i1)>SE(i2),
		ESI(k(c1)+1) = -1;
	else
		ESI(k(c1)) = -1;
		end;
	end;
if ~isempty(ESI),
	ESI  = ESI(ESI~=-1);	
	end;

NES = length(ESI);                 % Number of Energy Spikes

fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);

if NES<NSM,
	fprintf('Breaking because of insufficient spikes detected.	\n');
	return;
	end;

%=========================================================
% Spike Energy
% Plot Histogram of Peaks
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	h = Histogram(EP,Step,Step);
	hold on;
		h = plot(ETh*[1 1],[0 1000],'k:');
		set(h,'LineWidth',2*LWM);
		for c1 = 1:size(B,1),
			h1 = plot(B(c1,1),BD(c1,1),'g.');
			h2 = plot(B(c1,2),BD(c1,2),'r.');
			h3 = plot(B(c1,3),BD(c1,3),'g.');
			h  = [h1;h2;h3];
			set(h,'MarkerSize',20*MSM);	
			end;
		hold off;
	xlabel('Spike Peak Energy (normalized)');
	ylabel('Estimated PDF');
	title([NAME ': Spike Energy   Peak Distribution & Threshold']);	
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'SEPeakDistribution');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;

%=========================================================
% Spike Energy
% Plot Sample of Detected Spikes
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	L1  = max(round(ESI(1)-0.05*FS),1);
	%L2  = ESI(1)+0.20*FS;
	L2  = min(round(ESI(min(20,length(ESI)))+0.05*FS),NX);
	t2  = (0:NX-1)/FS;
	k   =  L1:L2;
	t   = (k-1)/FS;
	tmn = (L1-1)/FS;
	tmx = (L2-1)/FS;
	h   = plot(t,X(k),'b');
	set(h,'Color',[0.3 0.3 0.8]);
	hold on;
		h = plot(t,SE(k),'r');
		set(h,'Color',[1.0 0.2 0.2]);
		set(h,'LineWidth',1.2*LWM);
		h = plot(t2(PI(PI>=L1 & PI<=L2)),SE(PI(PI>=L1 & PI<=L2)),'r.');
		set(h,'MarkerSize',12*MSM);
		h = plot(t2(ESI(ESI>=L1 & ESI<=L2)),SE(ESI(ESI>=L1 & ESI<=L2)),'go');
		set(h,'Color'          ,GREEN);
		set(h,'MarkerFaceColor',GREEN);
		set(h,'MarkerSize',5*MSM);
		h = plot([tmn tmx],[ETh ETh],'k:');
		set(h,'LineWidth',2*LWM);
		hold off;
	axis([tmn tmx 1.1*min(X(k)) 1.1*max(X(k))]);
	xlabel('Time (sec)');
	ylabel('Spike Energy & Normalized Signal');
	title([NAME ': Spike Energy   Sample of Detected Spikes']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'SEDetectedSpikes');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;

%=========================================================
% Spike Energy
% Plot Worst Detected & Best Undetected 
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	ns = 15;
	k  = -NEP:NEP;
	t  = 1000*k/FS; % Time Index (ms)

	[jnk,in] = sort(SE(ESI));
	in       = in(ESI(in)+min(k)>=1 & ESI(in)+max(k)<=NX);
	nr       = min(ns,length(in));
	nc       = length(k);
	di       = diag(ESI(in(1:nr)))*ones(nr,nc) + ones(nr,nc)*diag(k); % Detected Indices

	ui       = PI(SE(PI)<=ETh);
	[jnk,in] = sort(-SE(ui));	
	in       = in(ui(in)+min(k)>=1 & ui(in)+max(k)<=NX);	
	nr       = min(ns,length(in));
	ui       = diag(ui(in(1:nr)))*ones(nr,nc) + ones(nr,nc)*diag(k); % Undetected Indices
	
	h = plot(t,X(ui)','r'); % Undetected in red
	hold on;	
		h = plot(t,X(di),'g');
		set(h,'Color',GREEN); % Detected in green
		hold off;
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	st = sprintf('Spike Energy   %d Worst Detected (Green)  %d Best Undetected (Red)',size(di,1),size(ui,1));
	title([NAME ': ' st]);	
	tmn = min(t);
	tmx = max(t);
	ymn = min(min([X(di);X(ui)]));
	ymx = max(max([X(di);X(ui)]));
	yrg = ymx-ymn;
	axis([tmn tmx ymn-0.05*yrg ymx+0.05*yrg]);
	set(gca,'Box','off');
	axisset(AFS);	
	drawnow;
  if PRINT==1,	
		st = sprintf('print %s %s;',PT,'SESignalDetectionBoundary');
		eval(st);
	elseif PRINT==2,
		print; 
  	end;

	figure;
	FigureSet(1,FGWD,FGHT);
	h = plot(t,SE(ui)','r'); % Undetected in red
	hold on;	
		h = plot(t,SE(di),'g');
		set(h,'Color',GREEN); % Detected in green
		h = plot([min(t) max(t)],[ETh ETh],'k:');
		set(h,'LineWidth',2*LWM);
		hold off;
	xlabel('Time (ms)');
	ylabel('Normalized Energy');
	st = sprintf('Spike Energy   %d Worst Detected (Green) and %d Best Undetected (Red)',size(di,1),size(ui,1));
	title([NAME ': ' st]);	
	tmn = min(t);
	tmx = max(t);
	ymn = min(min([SE(di);SE(ui)]));
	ymx = max(max([SE(di);SE(ui)]));
	yrg = ymx-ymn;
	axis([tmn tmx ymn-0.05*yrg ymx+0.05*yrg]);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,	
		st = sprintf('print %s %s;',PT,'SEDetectionBoundary');
		eval(st);
	elseif PRINT==2,
		print; 
  	end;
	end;

%=========================================================
% Spike Energy
% Plot the Overlapping Spikes
%========================================================= 
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	k  = -NEP:NEP;
	t = 1000*k/FS;
	h = [];
	si = ESI(ESI+min(k)>0 & ESI+max(k)<NX);
	hold on;
	ymn = mean(X);
	ymx = ymn;
	for cnt = 1:length(si),
   	hp  = plot(t,X(si(cnt)+k),'g');
   	h   = [h;hp];
		ymn = min(ymn,min(X(si(cnt)+k)));
		ymx = max(ymx,max(X(si(cnt)+k)));
		end;
	hold off;
	set(h,'Color',GREEN);
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	title([NAME ': Spike Energy   Overlapping Detected Spikes']);
	xmn = min(t);
	xmx = max(t);
	yrg = ymx-ymn;
	ymn = ymn-0.05*yrg;
	ymx = ymx+0.05*yrg;
	axis([xmn xmx ymn ymx]);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
   	st = sprintf('print %s %s;',PT,'SEOverlappingSpikes');
   	eval(st);
   elseif PRINT==2,
      print; 
   	end;
   end;

%=========================================================
% Template
% Construct Template
%=========================================================
fprintf('Building spike template...');
t0 = clock;

TDL = TDL/1000;                  % Convert from ms to sec
TOI = -TDL/2:TDL/(NTP-1):TDL/2;  % Calculate nominal locations of indices
TOI = TOI + OTP*(TDL/2);         % Offset template 
TOI = round(TOI*FS);             % Convert points to integers
TOI = unique(TOI);               % Remove redundant indices
NTP = length(TOI);               % Number of Template Points

TMP = zeros(size(TOI));
TP  = zeros(NES,NTP);

if NES==1,
	TP  = X(ESI+TOI);
	TMP = TP;
else
	% Phase I
	si = ESI;
	si = si(si+min(TOI)>=1 & si+max(TOI)<=NX);
	ni = length(si);                   % Number of indices
	for cnt = 1:ni,
		TP(cnt,:) = X(si(cnt)+TOI);
		end;
  	TMP = median(TP(1:ni,:));
    SETP  = TP;   
    SETMP = TMP;  

	% Phase II
	for c1 = 1:NTL-1,
		NSE      = TemplateError(X,TMP,TOI);
		mi       = FindPeaks(-NSE); % Minimum Indices
		[jnk,in] = sort(NSE(mi));
		np       = min(NES,NTS);
		mi       = mi(in(1:np));  % Only consider NTS or NES of them (NES = No. Energy Spikes)
		mi       = mi(NSE(mi)<median(NSE)); % Only use peaks that are less than median

		mi = mi(mi+min(TOI)>=1 & mi+max(TOI)<=NX);
		ni = length(mi); % Number of indices
		for cnt = 1:ni,
			TP(cnt,:) = X(mi(cnt)+TOI);
			end;

		TMP = median(TP(1:ni,:));
		end;
	TP = TP(1:ni,:);
	end;

fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);

%=========================================================
% Spike Energy Template
% Plot the Spike Energy Template and Spikes
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	t = 1000*TOI/FS;

	h = plot(t,SETP,'g');
	set(h,'Color',GREEN);
	hold on;
		h = plot(t,SETMP','b');
	 	set(h,'LineWidth',4*LWM);
		set(h,'Color',BROWN);
   	hold off;
	xmn = min(t);
	xmx = max(t);
	yrg = max(SETMP)-min(SETMP);
	ymn = min(SETMP)-0.2*yrg;
	ymx = max(SETMP)+0.2*yrg;
	axis([xmn xmx ymn ymx]);
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	title([NAME ': Spike Energy Template']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
   	st = sprintf('print %s %s;',PT,'SETemplateSpikes');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;

%=========================================================
% Template
% Plot the Template and Spikes Composing Template
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	t = 1000*TOI/FS;

	h = plot(t,TP,'g');
	set(h,'Color',GREEN);
	hold on;
		h = plot(t,TMP','b');
	 	set(h,'LineWidth',4*LWM);
		set(h,'Color',BROWN);
   	hold off;
	xmn = min(t);
	xmx = max(t);
	yrg = max(TMP)-min(TMP);
	ymn = min(TMP)-0.2*yrg;
	ymx = max(TMP)+0.2*yrg;
	axis([xmn xmx ymn ymx]);
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	title([NAME ': Final Template']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
   	st = sprintf('print %s %s;',PT,'TemplateSpikes');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;
  
%=========================================================
% Template Error
% Calculate Error, Find minimums, and Pick Threshold
%=========================================================
if PLOT==1,
	fprintf('Finding the template error spikes...');
	t0 = clock;

	NSE  = TemplateError(X,TMP,TOI);

	TEMI = FindPeaks(-NSE);
	TEMI = TEMI(NSE(TEMI)<median(NSE)); % Only use peaks that are less than median

	TEP = NSE(TEMI);          % Template Peaks

	[B,BD,Step] = FindBumps(TEP,ATH,NPP,MTH);
	nb   = size(B,1);
	if nb>=2,
		TETh = B(1,3);            % Set Template Threshold at right end of 1st bump
	else
		TETh = B(1,1);
		end;
	TESI = TEMI(TEP<TETh);    % Template Spike Indices

	% Eliminate spikes that are too close together (<TSM sec apart)
	sm  = (TSM/1000)*FS;
	dsi = diff(TESI);
	k   = 1:length(dsi);
	k   = k(dsi<sm);
	for c1 = 1:length(k),
		i1 = TESI(k(c1)  );
		i2 = TESI(k(c1)+1);
		if NSE(i1)<NSE(i2),
			TESI(k(c1)+1) = -1;
		else
			TESI(k(c1)  ) = -1;
			end;
		end;
	if ~isempty(TESI),
		TESI  = TESI(TESI~=-1);	
		end;

	NTES = length(TESI);      % Number of Template Spikes

	fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);
    end;

%=========================================================
% Template Error
% Plot Histogram of Minimums
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	h = Histogram(TEP,Step,Step);
	hold on;
		h = plot(TETh*[1 1],[0 1000],'k:');
		set(h,'LineWidth',2*LWM);
		for c1 = 1:size(B,1),
			h1 = plot(B(c1,1),BD(c1,1),'g.');
			h2 = plot(B(c1,2),BD(c1,2),'r.');
			h3 = plot(B(c1,3),BD(c1,3),'g.');
			h  = [h1;h2;h3];
			set(h,'MarkerSize',20*MSM);	
			end;
		hold off;
	xlabel('Template Error Minima (NSE)');
	ylabel('Estimated PDF');
	title([NAME ': Template Error   Peak Distribution & Threshold']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'TEDistribution');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;

%=========================================================
% Template Error
% Plot Sample of Detected Spikes
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	L1  = max(round(TESI(1)-0.05*FS),1);
	L2  = min(round(TESI(min(20,length(TESI)))+0.05*FS),NX);
	t2  = (0:NX-1)/FS;
	k   =  L1:L2;
	t   = (k-1)/FS;
	tmn = (L1-1)/FS;
	tmx = (L2-1)/FS;
	h   = plot(t,X(k),'b');
	set(h,'Color',[0.3 0.3 0.8]);
	hold on;
		h = plot(t,NSE(k),'r');
		set(h,'Color',[1.0 0.2 0.2]);
		set(h,'LineWidth',1.2*LWM);
		h = plot(t2(TEMI(TEMI>=L1 & TEMI<=L2)),NSE(TEMI(TEMI>=L1 & TEMI<=L2)),'r.');
		set(h,'MarkerSize',12*MSM);
		h = plot(t2(TESI(TESI>=L1 & TESI<=L2)),NSE(TESI(TESI>=L1 & TESI<=L2)),'go');
		set(h,'Color'          ,[0 1 0]);
		set(h,'MarkerFaceColor',[0 1 0]);
		set(h,'MarkerSize',5*MSM);
		h = plot([tmn tmx],[TETh TETh],'k:');
		set(h,'LineWidth',2*LWM);
		hold off;
	axis([tmn tmx 1.1*min(X(k)) median(NSE(k))]);
	xlabel('Time (sec)');
	ylabel('Template Error & Normalized Signal');
	title([NAME ': Template Error   Sample of Detected Spikes']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'TEDetectedSpikes');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;

%=========================================================
% Template Error
% Plot Worst Detected & Closest Undetected
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	ns = 15;
	k  = -NTP:NTP;
	t  = 1000*k/FS; % Time Index (ms)

	[jnk,in] = sort(-NSE(TESI));
	in       = in(TESI(in)+min(k)>=1 & TESI(in)+max(k)<=NX);
	nr       = min(ns,length(in));
	nc       = length(k);
	di       = diag(TESI(in(1:nr)))*ones(nr,nc) + ones(nr,nc)*diag(k); % Detected Indices

	ui       = TEMI(NSE(TEMI)>=TETh);
	[jnk,in] = sort(NSE(ui));	
	in       = in(ui(in)+min(k)>=1 & ui(in)+max(k)<=NX);	
	nr       = min(ns,length(in));
	ui       = diag(ui(in(1:nr)))*ones(nr,nc) + ones(nr,nc)*diag(k); % Undetected Indices

	h = plot(t,X(ui)','r'); % Undetected in red
	hold on;	
		h = plot(t,X(di),'g');
		set(h,'Color',GREEN); % Detected in green
		h = plot(1000*TOI/FS,TMP,'b');
		set(h,'LineWidth',4*LWM);
		set(h,'Color',BROWN);
		hold off;
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	st = sprintf('Template Error   %d NSE Worst Detected (Green)  %d Best Undetected (Red)',size(di,1),size(ui,1));
	title([NAME ': ' st]);	
	tmn = min(t);
	tmx = max(t);
	ymn = min([min([X(di);X(ui)]) TMP]);
	ymx = max([max([X(di);X(ui)]) TMP]);
	yrg = ymx-ymn;
	axis([tmn tmx ymn-0.05*yrg ymx+0.05*yrg]);
	set(gca,'Box','off');
	axisset(AFS);	
	drawnow;
  if PRINT==1,	
		st = sprintf('print %s %s;',PT,'TESignalDetectionBoundary');
		eval(st);
	elseif PRINT==2,
		print; 
  	end;

	figure;
	FigureSet(1,FGWD,FGHT);

	h = plot(t,NSE(ui)','r'); % Undetected in red
	hold on;	
		h = plot(t,NSE(di),'g');
		set(h,'Color',GREEN); % Detected in green
		h = plot([min(t) max(t)],[TETh TETh],'k:');
		set(h,'LineWidth',2*LWM);
		hold off;
	xlabel('Time (ms)');
	ylabel('Template Error (NSE)');
	st = sprintf('Template Error   %d Worst Detected (Green)  %d Best Undetected (Red)',size(di,1),size(ui,1));
	title([NAME ': ' st]);	
	tmn = min(t);
	tmx = max(t);
	ymn = min(min([NSE(di);NSE(ui)]));
	ymx = max(max([NSE(di);NSE(ui)]));
	yrg = ymx-ymn;
	axis([tmn tmx ymn-0.05*yrg ymx+0.05*yrg]);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,	
		st = sprintf('print %s %s;',PT,'TEDetectionBoundary');
		eval(st);
	elseif PRINT==2,
		print; 
  	end;
	end;

%=========================================================
% Template Error
% Plot the Overlapping Spikes
%========================================================= 
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	k  = -NTP:NTP;
	t  = 1000*k/FS;
	h  = [];
	si = TESI(TESI-min(k)>1 & TESI+max(k)<NX);
 	hold on;  
	ymn = mean(X);
	ymx = ymn;
	for cnt = 1:length(si),
   	hp  = plot(t,X(si(cnt)+k),'g');
   	h   = [h;hp];
		ymn = min(ymn,min(X(si(cnt)+k)));
		ymx = max(ymx,max(X(si(cnt)+k)));
		end;
	set(h,'Color',GREEN);
	h = plot(1000*TOI/FS,TMP,'b');
	set(h,'LineWidth',4*LWM);
	set(h,'Color',BROWN);
	hold off;
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	title([NAME ': Template Error   Overlapping Detected Spikes']);
	xmn = min(t);
	xmx = max(t);
	yrg = ymx-ymn;
	ymn = ymn-0.05*yrg;
	ymx = ymx+0.05*yrg;
	axis([xmn xmx ymn ymx]);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
   	st = sprintf('print %s %s;',PT,'TEOverlappingSpikes');
   	eval(st);
   elseif PRINT==2,
      print; 
   	end;
   end;

%=========================================================
% Template Correlation
% Calculate Correlation, Find Peaks & Pick Threshold
%=========================================================
fprintf('Finding the template correlation spikes...');
t0 = clock;

COR  = TemplateCorrelation(X,TMP,TOI);

TCPI = FindPeaks(COR);

TCP = COR(TCPI);        % Template Correlation Peaks

[B,BD,Step] = FindBumps(TCP,ATH,NPP,MTH);
nb    = size(B,1);      % Number of bumps 
if nb>=2,
	TCTh  = B(nb-1,3);      % Set Template Correlation Threshold
else
	TCTh  = B(1,3);  
	end;
TCSI  = TCPI(TCP>TCTh); % Template Correlation Spike Indices

% Eliminate spikes that are too close together (<TSM sec apart)
sm  = (TSM/1000)*FS;
dsi = diff(TCSI);
k   = 1:length(dsi);
k   = k(dsi<sm);
for c1 = 1:length(k),
	i1 = TCSI(k(c1)  );
	i2 = TCSI(k(c1)+1);
	if COR(i1)>COR(i2),
		TCSI(k(c1)+1) = -1;
	else
		TCSI(k(c1)) = -1;
		end;
	end;
if ~isempty(TCSI),
	TCSI  = TCSI(TCSI~=-1);	
	end;

NTCS  = length(TCSI);   % Number of Template Correlation Spikes

fprintf('done. (%5.2f min)\n',etime(clock,t0)/60);

%=========================================================
% Template Correlation
% Histogram of Correlation Peaks
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	h = Histogram(TCP,Step,Step);
	hold on;
		h = plot(TCTh*[1 1],[0 1000],'k:');
		set(h,'LineWidth',2*LWM);
		for c1 = 1:size(B,1),
			h1 = plot(B(c1,1),BD(c1,1),'g.');
			h2 = plot(B(c1,2),BD(c1,2),'r.');
			h3 = plot(B(c1,3),BD(c1,3),'g.');
			h  = [h1;h2;h3];
			set(h,'MarkerSize',20*MSM);	
			end;
		hold off;
	xlabel('Template Correlation');
	ylabel('Estimated PDF');
	title([NAME ': Template Correlation   Peak Distribution & Threshold']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'TCDistribution');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;

%=========================================================
% Template Correlation
% Plot Sample of Detected Spikes
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	L1  = max(round(TCSI(1)-0.05*FS),1);
	L2  = min(round(TCSI(min(20,length(TCSI)))+0.05*FS),NX);
	t2  = (0:NX-1)/FS;
	k   =  L1:L2;
	t   = (k-1)/FS;
	tmn = (L1-1)/FS;
	tmx = (L2-1)/FS;
	x   = X(k)*2/max(abs(X(k)));
	h   = plot(t,x,'b');
	set(h,'Color',[0.3 0.3 0.8]);
	hold on;
		h = plot(t,COR(k),'r');
		set(h,'Color',[1.0 0.2 0.2]);
		set(h,'LineWidth',1.2*LWM);
		h = plot(t2(TCPI(TCPI>=L1 & TCPI<=L2)),COR(TCPI(TCPI>=L1 & TCPI<=L2)),'r.');
		set(h,'MarkerSize',12*MSM);
		h = plot(t2(TCSI(TCSI>=L1 & TCSI<=L2)),COR(TCSI(TCSI>=L1 & TCSI<=L2)),'go');
		set(h,'Color'          ,GREEN);
		set(h,'MarkerFaceColor',GREEN);
		set(h,'MarkerSize',5*MSM);
		h = plot([tmn tmx],[TCTh TCTh],'k:');
		set(h,'LineWidth',2*LWM);
		hold off;
	axis([tmn tmx 1.1*min(x) 1.1*max(x)]);
	xlabel('Time (sec)');
	ylabel('Template Correlation & Normalized Signal (Scaled)');
	title([NAME ': Template Correlation   Sample of Detected Spikes ']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'TCDetectedSpikes');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;

%=========================================================
% Template Correlation
% Plot Worst Detected & Closest Undetected
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	ns = 15;
	k  = -NTP:NTP;
	t  = 1000*k/FS; % Time Index (ms)

	[jnk,in] = sort(COR(TCSI));
	in       = in(TCSI(in)+min(k)>=1 & TCSI(in)+max(k)<=NX);
	nr       = min(ns,length(in));
	nc       = length(k);
	di       = diag(TCSI(in(1:nr)))*ones(nr,nc) + ones(nr,nc)*diag(k); % Detected Indices

	ui       = TCPI(COR(TCPI)<=TCTh);
	[jnk,in] = sort(-COR(ui));
	in       = in(ui(in)+min(k)>=1 & ui(in)+max(k)<=NX);
	nr       = min(ns,length(in));
	ui       = diag(ui(in(1:nr)))*ones(nr,nc) + ones(nr,nc)*diag(k); % Undetected Indices

	h = plot(t,X(ui)','r'); % Undetected in red
	hold on;	
		h = plot(t,X(di),'g');
		set(h,'Color',GREEN); % Detected in green
		h = plot(1000*TOI/FS,TMP,'b');
		set(h,'LineWidth',4*LWM);
		set(h,'Color',BROWN);
		hold off;
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	st = sprintf('Template Correlation   %d Worst Detected (Green)  %d Best Undetected (Red)',size(di,1),size(ui,1));
	title([NAME ': ' st]);	
	tmn = min(t);
	tmx = max(t);
	ymn = min([min([X(di);X(ui)]) TMP]);
	ymx = max([max([X(di);X(ui)]) TMP]);
	yrg = ymx-ymn;
	axis([tmn tmx ymn-0.05*yrg ymx+0.05*yrg]);
	set(gca,'Box','off');
	axisset(AFS);	
	drawnow;
  if PRINT==1,	
		st = sprintf('print %s %s;',PT,'TCSignalDetectionBoundary');
		eval(st);
	elseif PRINT==2,
		print; 
  	end;

	figure;
	FigureSet(1,FGWD,FGHT);

	h = plot(t,COR(ui)','r'); % Undetected in red
	hold on;	
		h = plot(t,COR(di),'g');
		set(h,'Color',GREEN); % Detected in green
		h = plot([min(t) max(t)],[TCTh TCTh],'k:');
		set(h,'LineWidth',2*LWM);
		hold off;
	xlabel('Time (ms)');
	ylabel('Template Correlation');
	st = sprintf('Template Correlation   %d Worst Detected (Green)  %d Best Undetected (Red)',size(di,1),size(ui,1));
	title([NAME ': ' st]);	
	tmn = min(t);
	tmx = max(t);
	ymn = min(min([COR(di);COR(ui)]));
	ymx = max(max([COR(di);COR(ui)]));
	yrg = ymx-ymn;
	axis([tmn tmx ymn-0.05*yrg ymx+0.05*yrg]);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,	
		st = sprintf('print %s %s;',PT,'TCDetectionBoundary');
		eval(st);
	elseif PRINT==2,
		print; 
  	end;
	end;

%=========================================================
% Template Correlation
% Plot the Overlapping Spikes
%========================================================= 
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	k  = -NTP:NTP;
	t  = 1000*k/FS;
	h  = [];
	si = TCSI(TCSI-min(k)>1 & TCSI+max(k)<NX);
 	hold on;  
	ymn = mean(X);
	ymx = ymn;
	for cnt = 1:length(si),
   	hp  = plot(t,X(si(cnt)+k),'g');
   	h   = [h;hp];
		ymn = min(ymn,min(X(si(cnt)+k)));
		ymx = max(ymx,max(X(si(cnt)+k)));
		end;
	set(h,'Color',GREEN);
	h = plot(1000*TOI/FS,TMP,'b');
	set(h,'LineWidth',4*LWM);
	set(h,'Color',BROWN);
	hold off;
	xlabel('Time (ms)');
	ylabel('Normalized Signal');
	title([NAME ': Template Correlation   Overlapping Detected Spikes']);
	xmn = min(t);
	xmx = max(t);
	yrg = ymx-ymn;
	ymn = ymn-0.05*yrg;
	ymx = ymx+0.05*yrg;
	axis([xmn xmx ymn ymx]);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
   	st = sprintf('print %s %s;',PT,'TCOverlappingSpikes');
   	eval(st);
   elseif PRINT==2,
      print; 
   	end;
   end;

%=========================================================
% Resolve Differences between Energy Spikes & Error Spikes
%========================================================= 
%SI = TESI;
%SI = ESI;
SI = TCSI;

if length(SI)<NSM,
	fprintf('Breaking because of insufficient spikes detected.	\n');
	return;
	end;

ST = (TCSI-1)/FS;

%=========================================================
% Plot a Sample of the Spikes
%=========================================================
if PLOT==1,
	figure;
	FigureSet(1,FGWD,FGHT);

	L1  = max(round(SI(1)-0.05*FS),1);
	%L2  = SI(1)+0.20*FS;
	L2  = min(round(SI(min(20,length(SI)))+0.05*FS),NX);
	k   =  L1:L2;
	t   = (0:NX-1)/FS;
	ymax = 1.1*max(X(k));
	ymin = 1.1*min(X(k)); 
	si   = SI(SI<=max(k));
	y    = ymax*ones(size(si));
	h1   = stem(t(si),y,'r:');
	hold on;
		y  = ymin*ones(size(si));
		h2 = stem(t(si),y,'r:');   
		h  = [h1;h2];  
		set(h,'Marker','none');
		set(h,'LineWidth',1.5*LWM);
		plot(t(k),X(k),'b');
		hold off;
	axis([min(t(k)) max(t(k)) ymin ymax]);
	xlabel('Time (sec)');
	ylabel('Normalized Signal');
	title([NAME ': Sample of Detected Spikes']);
	set(gca,'Box','off');
	axisset(AFS);
	drawnow;
	if PRINT==1,
		st = sprintf('print %s %s;',PT,'DetectedSample');
		eval(st);
	elseif PRINT==2,
		print; 
		end;
	end;
   
