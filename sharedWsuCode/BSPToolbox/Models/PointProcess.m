function [b,ei] = PointProcess(x,fsa,fra,mda,sha,tra,pfa);
%PointProcess: Creates a binary event series from a given intensity
%
%   [b,si] = PointProcess(x,fs,fr,md,sh,tr,pf);
%
%   x    Process intensity. Scaled (-1<x<1).
%   fs   Sample rate (Hz). Default = 1 Hz.
%   fr   Mean firing rate. Default = 50 Hz.
%   md   Modulation index. Default = 0.8.
%   sh   Gamma function shape paramter. Default = 1.
%   tr   Refractory time. Default = 0.001 s.
%   pf   Plot flag: 0=none (default), 1=screen, 2=current figure.
%
%   b    Binary event series.
%   ei   Intervals between events.
%
%   Creates a synthetic binary event series with the specified
%   intensity using the model given in the citation below.
%
%   Example: Generate a binary event series with a tremor intensity.
%
%        fm = 6;                             % Mean tremor frequency
%        fs = 1000;                          % Sample rate (Hz)
%        ny = 60*fs;                         % Duration
%        wn = decimate(10*randn(ny+1e3,1),50);
%        fi = Lowpass(wn,fs/50,0.5);
%        fi = resample(fi,50,1);
%        fi = fi(end-ny+1:end) + fm;         % Eliminate edge effect
%        th = cumsum(fi*2*pi/fs);            % Instantaneous angle
%        it = sin(th);                       % Quasi-periodic sinusoid
%        b = PointProcess(it,fs,50,0.9,1);
%        NonparametricSpectrogram(decimate(b,20),fs/20,1);
%        hold on; plot(((1:ny)-0.5)/fs,fi,'w'); hold off;
%
%   J. McNames, "Pulse Frequency Demodulation Without Spike
%   Detection," Proceedings of the 2nd International IEEE EMBS
%   Conference on Neural Engineering, pp. 229-232, March 2005.
%
%   Version 1.00 JM
%
%   See also KernelFilter and KalmanFilterFrequencyTracker.

%====================================================================
% Error Checking
%====================================================================
if nargin<2,
   help PointProcess;
   return;
   end;

if max(x)>1 | min(x)<-1,
   error('Requires intensity signal input bounded between -1 and 1.');
   end;

%=====================================================================
% Process function arguments
%=====================================================================
fs = 1;
if exist('fsa','var') & ~isempty(fsa),
   fs = fsa;
   end;

fr = 50;                                                   % Mean firing rate
if exist('fra','var') & ~isempty(fra),
   fr = fra;
   end;

tr = 0.001;                                                % Refractory period
if exist('tra','var') & ~isempty(tra),
   tr = tra;
   end;

md = 0.8;                                                  % Modulation index
if exist('mda','var') & ~isempty(mda),
   md = mda;
   end;

sh = 1.0;                                                  % Shape parameter
if exist('sha','var') & ~isempty(sha),
   sh = sha;
   end;

pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
   pf = 1;
   end;
if exist('pfa') & ~isempty(pfa),
   pf = pfa;
   end;

%====================================================================
% Pre-processing
%====================================================================
x   = x(:);                                                % Make into a column vector
nx  = length(x);                                           % Duration (samples)
tx  = length(x)/fs;                                        % Duration (seconds)
nr  = round(tr*fs);                                        % Refractory period (samples)

%====================================================================
% Generate the Signal Intensity
%====================================================================
it = fr + fr*md*x;                                         % Signal intensity (spikes/second) = mean firing rate + quasi-periodic firing rate
sit = it/(1-fr*tr);                                        % Scale the intensity to account for refractory effects

%====================================================================
% Generate the Intervals
%====================================================================
cit = cumsum(sit)/fs;                                      % Calculate the piecewise intervals
b   = zeros(nx,1);                                         % Spike train
ei  = zeros(ceil(10*tx/fr),1);                             % Allocate memory for the spike intervals
ni  = 0;                                                   % Initialize a counter
i0  = 1;                                                   % Initial index
i1  = i0+nr;                                               % Index after a refractory time
while sum(ei(1:ni))<tx,
   e = gamrnd(sh,1/sh);
   while i1<=nx && i0+nr<=nx && cit(i1)-cit(i0+nr)<e,     % Core process to determine whether or not
       i1 = i1 + 1;                                       % To generate action potentials
       end;
   if i1>nx || i0+nr>nx,
       break;
       end;
   ni      = ni + 1;
   ei(ni)  = (i1-i0)/fs;
   b(i1-1) = 1;
   i0      = i1;
   i1      = i0+nr;
   %disp([e sum(ei) tx]); drawnow;
   end;
ei = ei(1:ni);

%==================================================================
% Plot the Results
%==================================================================
if pf>=1,
   figure;
   FigureSet(2);
   k = 1:nx;
   t = (k-0.5)/fs;
   h = stem(t,b);
   set(h,'Marker','None');
   hold on;
       h = plot(t,it/fs,'r');
       hold off;
   xlabel('Time (s)');
   ylabel('Binary Event Series');
   box off;
   AxisSet;

   ds = floor(fs/(2*fr));
   NonparametricSpectrogram(decimate(b,ds),fs/ds,30/fr,[],[],[],pf);
   end;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
   clear('b','ei');
   end;
