function [M] = DensityModes(x,para,mpa,mda,pfa)
%DensityModes: Identifies dominant modes in multimodal distributions. 
%
%   [M] = DensityModes(x,ar,mp,md,pf);
%
%   x     Vector containing the univariate data (unsorted).
%   par   Required area for a peak to be considered a mode. Default = 5e-3. 
%   mp    Maximum number of modes to locate. Default = 4.
%   md    Minimum depth required between two density peaks for the
%         peaks to qualify as modes (fraction of smallest peak height).
%         Default = 0.75.
%
%   M     Matrix containing mode characteristics, 1 row per mode.
%   fc    Cutoff frequency of the 
%  
%   Uses density estimation to locate modes in the data. The algorithm
%   initially oversmoothes the data and then reduces the smoothness 
%   until the distinct modes are identified.
%
%   The matrix of mode characteristics contains seven columns. The 
%   contents of each column is as follows:
%
%      1: Location of the preceding (left of the mode) density minimum
%      2: Location of the density peak
%      3: Location of the following (right of the mode) density minimum
%      4: Area of the mode, as segmented by the minima
%      5: Density at the preceding minimum
%      6: Density at the peak 
%      7: Density at the following minimum
%
%   Example: Find the modes in the marginal density of a noisy 
%   electrocardiogram.
%
%      load Tremor;
%      [pi,y] = PowerPeaks(x,fs,200);
%      [M] = DensityModes(y(pi),[],[],[],1);
%
%   Desmond J. Higham and Nicholas J. Higham, "MATLAB Guide," 2000.
%
%   Version 2.00.22 JM
%
%   See also DetectMaxima, Lowpass, and ECGDetectREnergy.

%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1,
    help DensityModes;
    return;
    end;

par = 5e-3;                            % Required peak area
if exist('para') & ~isempty(para),
    par = para;
    end;
 
mp = 4;                                % Maximum no. peaks
if exist('mpa') & ~isempty(mpa),
    mp = mpa;
    end;
    
md = 0.75;                             % Minimum depth between peaks (fraction of smallest peak height)
if exist('mda') & ~isempty(mda),
    md = mda;
    end;

pf = 0;                                % Default - no plotting
if nargout==0,   % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
 
%=========================================================================
% Author-Specified Parameters
%=========================================================================
nb  = 50;   % No. of bins to use in first density estimate
nl  = 15;   % Maximum no. of loops (bin width reductions)
bwr = 0.80; % Rate of bin width reduction

%=========================================================================
% Preprocessing
%=========================================================================
x = x(:); % Make into a column vector
xmin = min(x);
xmax = max(x);

bw = (xmax-xmin)/nb; % Initial bin size

cm = 0; % Count of modes
M  = zeros(mp,7);

%=========================================================================
% Main Loop
%=========================================================================
for c1 = 1:nl,
    kw = bw*10;
	bc = ((xmin-5*kw):bw:(xmax+5*kw)).';   
    nb = length(bc);                         % No. bins
    yh = hist(x,bc);                         % Put the data in sparse, sorted bins
    yh = yh/(sum(yh)*bw);                    % Normalize to have unit area
    d  = SmoothSeries(bc,yh,bc,kw);          % Smooth the histogram to generate estimated density
	    
    k  = 2:nb-1;
    pi = k((d(k)>d(k-1) & d(k)>=d(k+1)) | (d(k)>=d(k-1) & d(k)>d(k+1))); % Peak indices  
    di = k((d(k)<d(k-1) & d(k)<=d(k+1)) | (d(k)<=d(k-1) & d(k)<d(k+1))); % Depression indices  
    
    if d(1)<=d(2),
        di = [1,di];
        end;
    if d(nb)<=d(nb-1),
        di = [di,nb];
        end;
    
    %---------------------------------------------------------------------
    % Check Maximum No. Peaks Criterion
    %---------------------------------------------------------------------    
    if c1~=1 & length(pi)>mp, % If more than the maximum no. of peaks is found, break
        break;
        end;    
        
    %---------------------------------------------------------------------
    % Identify peaks with sufficiently deep minima between them    
    %---------------------------------------------------------------------
    ci = 1; % Current index
    ri = 2; % Right index    
    while ri<=length(pi),
        while ri<=length(pi),
            mi = di(di>pi(ci) & di<pi(ri));  % Indices of minima between peaks
            if length(mi)>1,                 % If more than one, pick the smallest minimum
                [jnk,id] = min(d(mi));       
                mi = mi(id);
                end;
                
            if d(mi)>md*min([d(pi(ci)) d(pi(ri))]), % If minimum is insufficiently deep,
                if d(pi(ci))<d(pi(ri)),      %   If current peak is the smaller,
                    pi(ci) = -1;             %     Label peak as invalid
                    ci     = ri;             %     Move on to next valid peak
                    ri     = ri + 1;         %     Increment right peak index
                else                         %   Otherwise
                    pi(ri) = -1;             %     Label right peak as invalid
                    ri     = ri + 1;         %     Increment to next right peak
                    end; 
            else                             % Otherwise (minimum is sufficiently deep)
                ci = ri;                     %   Move on to next peak
                ri = ri + 1;                 %   Increment right peak index
                break;
                end;
            end;
        end;
    pi = pi(pi~=-1);  % Strip out the invalid peaks
    
    %---------------------------------------------------------------------
    % Identify peaks with sufficient area
    %---------------------------------------------------------------------
    cm = 0; % Count of modes
    for c2 = 1:length(pi), 
        % Find the left minimum
        if c2==1,
            id = 1;
        else
            id = pi(c2-1);
            end;
        li = di(di>=id & di<pi(c2));         % Index of left minimum (left edge of mode)
        if length(li)>1,                     % If more than one, pick the smallest minimum
            [jnk,id] = min(d(li));       
            li = li(id);
            end;
        
        % Find the right minimum
        if c2==length(pi),
            id = length(d);
        else
            id = pi(c2+1);
            end;
        ri = di(di> pi(c2) & di<=id);        % Index of right minimum (right edge of mode)
        if length(ri)>1,                     % If more than one, pick the smallest minimum
            [jnk,id] = min(d(ri));       
            ri = ri(id);
            end;        
        
        pa = sum(d(li:ri))*bw;
        
        if pa>par,                % If mode area is sufficient, save as a valid mode
            cm = cm + 1;
            M(cm,1) = bc(li);                % Location of left edge of mode
            M(cm,2) = bc(pi(c2));            % Location of peak of mode
            M(cm,3) = bc(ri);                % Location of right edge of mode
            M(cm,4) = pa;                    % Estimated area of mode
            M(cm,5) = d(li);                 % Estimated density at left edge of mode
            M(cm,6) = d(pi(c2));             % Estimated density at peak
            M(cm,7) = d(ri);                 % Estimated density at right edge of mode
            end;
        end;
        
    if 0,
        FigureSet;
        h = plot(bc,d);
        hold on;
        for c2 = 1:cm,
            h = plot(M(c2,1),M(c2,5),'g.');
            h = plot(M(c2,2),M(c2,6),'r.');
            h = plot(M(c2,3),M(c2,7),'g.');
            end;
            hold off;
        title('Estimated Probability Density Function (PDF)');
        xlabel('Amplitude');
        ylabel('PDF (scaled to unit area)');
        xlim([min(bc) max(bc)]);
        ylim([0 1.02*max(d)]);
        box off;
        zoom on;
        AxisSet;    
        fprintf('No. modes: %d. Pausing...\n',cm);
        pause;
        end;
        
    bw = bw*bwr;        
    
    bcs = bc; % Bin Centers - save for plotting
    ds  = d;  % Density values - save for plotting
    
    end;
    
M = M(1:cm,:);                               % Eliminate empty rows

%====================================================================
% Plot Density Plot with Bumps Labeled
%====================================================================
if pf==1,
    figure;
    FigureSet;
    h = plot(bcs,ds);
    set(h,'LineWidth',1.2);
    nm = size(M,1); % No. modes
    hold on;
    for c1 = 1:nm,
        h1 = plot(M(c1,1),M(c1,5),'g.');
        h2 = plot(M(c1,2),M(c1,6),'r.');
        h3 = plot(M(c1,3),M(c1,7),'g.');
        h  = [h1;h2;h3];
        set(h,'MarkerSize',15);
        end;
    hold off;
    title('Estimated Probability Density Function (PDF)');
    xlabel('Amplitude');
    ylabel('PDF (scaled to unit area)');
    xlim([min(bcs) max(bcs)]);
    ylim([0 1.02*max(d)]);
    box off;
    zoom on;
    AxisSet;
    end

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('M');
    end;

    