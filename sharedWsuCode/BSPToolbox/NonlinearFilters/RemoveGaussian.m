function [Y] = RemoveGaussian(X, wha, wva, dha, dva, pfa)
%RemoveGaussian: Two-Dimensional Gaussian Noise Removing Filter
%
%   [Y] = RemoveGaussian(X,wh,wv,dh,dv,pf)
%
%   X       Input signal(must be a matrix)
%   wh      Horizontal size of the sliding window in samples (odd 
%           integer). Default=21.
%   wv      Vertical size of the sliding window in samples (odd 
%           integer). Default=21.
%   dh      Size of the horizontal step of the sliding window (samples).
%           Default=1.
%   dv      Size of the vertical step of the sliding window (samples). 
%           Default=1.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   Y       Filtered signal
%
%   Filters a two-dimensional signal using a Mean Filter. The matrix
%   is padded by repeating the values at each edge. Then, a sliding
%   window of specified width and length is placed at the initial value
%   of the input matrix. The mean of the values inside the window is 
%   calculated and stored in the output matrix. The same procedure is 
%   repeated as the window slides through the data. The number of 
%   the window advances at each step in the horizontal and vertical 
%   directions is determined by the values of the input parameters dh
%   and dv, respectively.
%
%   Example: Filter the spectrogram of the ICP signal using a Gaussian
%   Noise Removing Filter with a window of 11-by-21 samples.
%
%       load ICP.mat;
%       icpd = decimate(icp, 15);
%       [S,t,f] = Spectrogram(icpd,125/15,[],[],[],[],1);
%       [Y] = RemoveGaussian(S, 11,21,5,5,1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, 1997.
%
%   Version 1.00 CC
%
%   See also RemoveImpulses. 

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>6,
    help RemoveGaussian;
    return;
    end;

wh = 21; % Default window horizontal size = 21 samples
if exist('wha') & ~isempty(wha),
    wh = wha;
    end;

wv = 21; % Default window vertical size = 21 samples
if exist('wva') & ~isempty(wva),
    wv = wva;
    end;

dh = 1; % Default window horizontal step = 1 sample
if exist('dha') & ~isempty(dha),
    dh = dha;
    end;

dv = 1; % Default window vertical size = 1 sample
if exist('dva') & ~isempty(dva),
    dv = dva;
    end;
    
pf = 0; % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

% -------------------------------------------------------------------
% Create padded matrix
% -------------------------------------------------------------------
sx = size(X);                   % Size of input matrix X
hx = sx(1);                     % Horizontal size (number of rows)
vx = sx(2);                     % Vertical size (number of columns)
ph = (wh-1)/2;                  % Horizontal Pad size
pv = (wv-1)/2;                  % Vertical Pad Size

Xp = zeros(hx+2*ph,vx+2*pv);    % Create padded matrix of zeros


for i=1:hx
    Xp(ph+i,:) = [X(i,1)*ones(1,pv) X(i,:) X(i,vx)*ones(1, pv)];
end

for j=1:vx+2*pv
    Xp(:,j) = [Xp(ph+1,j)*ones(ph,1) ; Xp(ph+1:ph+hx,j) ; Xp(hx+ph,j)*ones(ph,1)]; 
end

%--------------------------------------------------------------------
% Filter signal
%--------------------------------------------------------------------
sp = size(Xp);
hp = sp(1);
vp = sp(2);


ci = 1;
cj = 1;
for i = ph+1:dh:hp-ph
    for j = pv+1:dv:vp-pv
        Xw = Xp(i-ph:i+ph,j-pv:j+pv);
        sXw = size(Xw);
        Xw = reshape(Xw,sXw(1)*sXw(2),1);
        Y(ci,cj) = mean(Xw);
        cj = cj+1;
    end
    ci = ci+1;
    cj = 1;
end

sy = size(Y);
hy = sy(1);
vy = sy(2);

%--------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------
if exist('pf') & pf == 1,
    figure
    ha1 = axes;
    s = reshape(abs(Y),vy*hy,1);
    p = [0 prctile(s,98)];
    imagesc(1:hy,1:vy,abs(Y),p); 
    xlim([0 hy]);
    ylim([0 vy]);
    set(ha1,'YDir','normal');
    xlabel('Samples');
    ylabel('Samples');
    title('Gaussian Noise Removing Filter');
    AxisSet(10); 
    FigureSet(1);
end
