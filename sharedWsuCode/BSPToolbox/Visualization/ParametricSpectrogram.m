function [S,t,f] = ParametricSpectrogram(x,fsa,wla,foa,eta,wta,fra,nfa,nsa,pfa);
%ParametricSpectrogram: Generates the spectrogram of the signal
%
%   [S,t,f] = ParametricSpectrogram(x,fs,wl,fo,et,wt,fr,nf,ns,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1024 samples.
%        If a vector, specifies entire window.
%   fo   Model order. Default = 30. 
%   et   Estimator type: 0=Unbiased/no window, 
%        1=Modified covariance (default). 
%   wt   Window type: 0=Rectangular, 1=Blackman (default).
%   fr   Minimum and maximum frequencies to display (Hz).
%        Default = [0 fs/2].
%   nf   Number of frequencies to evaluate (vertical resolution). 
%        Default = max(128,round(wl/2)).
%   ns   Requested number of times (horizontal pixels) to evaluate 
%        Default = min(400,length(x)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   S    Matrix containing the image of the Spectrogram.
%   t    Times at which the spectrogram was evaluated (s).
%   f    Frequencies at which the spectrogram was evaluated (Hz).
%
%   Calculates estimates of the spectral content at the specified 
%   times using an autoregressive model. The mean of the signal is 
%   removed as a preprocessing step. The square root of the power 
%   spectral density is calculated and displayed. To limit 
%   computation, decimate the signal if necessary to make the upper 
%   frequency range approximately equal to half the sampling 
%   frequency. 
%   
%   If only the window length is specified, the blackman window is 
%   applied as a data window. The specified window length should be 
%   odd. If called with a window with an even number of elements, a 
%   zero is appended to make the window odd.
%
%   Example: Generate the parametric spectrogram of an intracranial 
%   pressure signal using a Blackman-Harris window that is 45 s in 
%   duration.
%
%      load ICP.mat; 
%      x = decimate(icp,5);
%      ParametricSpectrogram(x,fs/5,30,30);
%
%   D. G. Manolakis, V. K. Ingle, S. M. Kogon, "Statistical and 
%   Adaptive Signal Processing," McGraw-Hill, 2000.  
%
%   Version 1.00 JM
%
%   See also specgram, window, decimate, and Spectrogram.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help ParametricSpectrogram;
    return;
    end;
    
if length(x)==0,
    error('Signal is empty.');
    end;

if var(x)==0,
    error('Signal is constant.');
    end;
    
if exist('fra') & length(fra)==1,
    error('Frequency range must be a 2-element vector.');
    end;
    
%====================================================================
% Calculate Basic Signal Statistics
%====================================================================    
nx = length(x);
xr = x;
mx = mean(x);
sx = std(x);
    
%====================================================================
% Process Function Arguments
%====================================================================
fs = 1;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;    
   
wl = min(1024,nx);                                         % Default window length
if exist('wla') & ~isempty(wla),
    wl = max(3,round(wla*fs));
    if ~rem(wl,2),
        wl = wl + 1;                                   % Make odd
        end;
    end;

fo = 30;                                                   % Default filter oder
if exist('foa') & ~isempty(foa),
    fo = foa;
    end;       

et = 1;                                                    % Modified covariance estimation type
if exist('eta') & ~isempty(eta),
    et = eta;
    end;       
    
wt = 1;                                                    % Blackman data matrix window
if exist('wta') & ~isempty(wta),
    wt = wta;
    end;       
    
fmin = 0;                                                  % Lowest frequency to display
fmax = fs/2;                                               % Highest frequency to display
if exist('fra') & ~isempty(fra),
    fmin = max(fra(1),0);
    fmax = min(fra(2),fs/2);
    end;        
    
nf = max(128,round(wl/2));
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
    end;   
    
ns = min(400,nx);
if exist('nsa') & ~isempty(nsa),
    ns = min(nsa,nx);
    end;         
    
pf = 0;                                                    % Default - no plotting
if nargout==0,                       % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% More Error Checking
%====================================================================    
if fo>wl,
    error('Filter order is larger than the window length.');
    end;
    
%====================================================================
% Preprocessing
%====================================================================  
switch wt,
case 0, 
    wn = ones(wl-fo,1);
case 1, 
    wn = blackman(wl-fo+2);                                   
    wn = wn(2:end-1);                                      % Trim off the zeros
    end;
x  = x(:).';                                               % Make into a row vector
wn = wn(:);                                                % Make into a column vector
x  = x - mx;                                               % Remove mean

%====================================================================
% Variable Allocations & Initialization
%====================================================================
st   = (nx-1)/(ns-1);                                      % Step size (samples)
st   = max(st,1);
te   = 1:st:nx;                                            % Evaluation times (in samples)
nt   = length(te);                                         % No. evaluation times
nz   = nf*(fs/(fmax-fmin));                                % No. window points needed in FFT
nz   = 2^(ceil(log2(nz)));                                 % Convert to power of 2 for FFT
b0   = floor(nz*(fmin/fs))+1;                              % Lower frequency bin index 
b1   = ceil (nz*(fmax/fs))+1;                              % Upper frequency bin index
fi   = b0:b1;                                              % Frequency indices
f    = fs*(fi-1)/nz;                                       % Frequencies
nf   = length(fi);                                         % No. of frequencies that PSD is evaluated at 

%====================================================================
% Main loop
%====================================================================
switch et,
case 0,                                                    % Short/no window
    X = zeros(wl-fo,fo);
    y = zeros(wl-fo,1);
case 1,                                                    % Modified Covariance
    X = zeros(2*(wl-fo),fo);                               
    y = zeros(2*(wl-fo),1);                        
    end;

S  = zeros(nf,nt);                                         % Power Spectral Density
for c1 = 1:nt,
    ic  = te(c1);                                          % Sample index of segment center
    i0  = round(ic-wl/2);                                  % Index of first segment element               
    i1  = i0 + (wl-1);                                     % Index of last segment element
    k   = (i0:i1);                                         % Collection of indices
    a0  = max(-i0+1,0);                                    % Number of zeros to append at front
    a1  = max(i1-nx,0);                                    % Number of zeros to append at end
    k   = k(1+a0:wl-a1);                                   % Subset of valid indice    
    xs  = [x(1)*ones(1,a0) x(k) x(nx)*ones(1,a1)].';       % Data segment        
 
    %----------------------------------------------------------------
    % Create the Target Vector and Data Matrix (full windowing)
    %----------------------------------------------------------------
    y(1:wl-fo) = wn.*xs(fo+1:wl); 
    for c2=1:fo,
       X(1:wl-fo,c2) = wn.*xs(fo-(c2-1):wl-c2);
       end;
    if et==1,
        y(wl-fo+(1:wl-fo)) = wn.*xs(1:wl-fo); 
        for c2=1:fo,
           X(wl-fo+(1:wl-fo),c2) = wn.*xs(c2+1:wl-fo+c2);
           end;        
        end;
   	ah  = -pinv(X)*y;
    s2w = mean(abs(y + X*ah).^2); 
    ah  = [1;ah];
    
    %[ah,s2w] = armcov(y,mo);
    
% 	xs  = y.*wn;
% 	
%     np  = 2^nextpow2(2*wl-1);                              % Figure out how much to pad the signal
%     XS  = fft(xs,np);                                      % Calculate the FFT
%     ac  = real(ifft(abs(XS).^2)).'/wl;                     % Calculate the inverse FFT
%     R   = toeplitz(ac(1:mo),ac(1:mo));                     % Create autocorrelation matrix
%     d   = ac(2:mo+1);                                      % Create cross-correlation vector
%     a   = R\d; 
%     ah  = [1; -a]; 
%     s2w = ac(1:fo+1)'*ah;                                  % Estimated variance of the excitation process w(n). Also equal to sum(xs.^2)/wl-d'*a; 
    
    dft = s2w./(abs(fft(ah,nz)).^2);                       % Calculate FFT
    
    S(:,c1) = dft(fi).';    
    drawnow;
	end;
t = (te-1)/fs;

%====================================================================
% Postprocessing
%====================================================================  
x = xr;   
t = t(:); % Convert to column vector
f = f(:); % Convert to column vector    

%====================================================================
% Plot Results
%====================================================================  
if pf>=1,
    td = 1; % Time Divider
    if max(t)>2000,
        td = 60; % Use minutes
        end;
    
    Tmax = (nx-1)/fs;
    
    if pf~=2,
        figure;
        FigureSet(1);
        colormap(jet(256));
        end;       
        
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
    % Parametric Spectrogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.15 0.21 0.81 0.69]);
    s = reshape(abs(S),nf*nt,1);
    
    p = [0 prctile(s,98)];
    imagesc(tp,f,abs(S),p); %,[Smin Smax]);

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
    ha2 = axes('Position',[0.135 0.21 0.01 0.69]);    
    colorbar(ha2);    
    set(ha2,'Box','Off')
    set(ha2,'YTick',[]);
    
    %----------------------------------------------------------------
    % Power Spectral Density
    %----------------------------------------------------------------
    ha3 = axes('Position',[0.08 0.21 0.05 0.69]);
    psd = mean((abs(S).^2).').^(1/2);
    plot(psd,f,'r'); %,[Smin Smax]);
    ylim([fmin fmax]);
    xlim([0 max(psd)]);
    ylabel('Frequency (Hz)');
    set(gca,'XAxisLocation','Top');    
    
    %----------------------------------------------------------------
    % Signal
    %----------------------------------------------------------------      
    ha4 = axes('Position',[0.15 0.10 0.81 0.10]);

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
    
    axes(ha1);
    AxisSet;
    end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('S','t','f');
    end;

