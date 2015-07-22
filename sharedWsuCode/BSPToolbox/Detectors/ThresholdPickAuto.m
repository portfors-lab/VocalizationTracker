function [th,thmax] = PickThreshold(x,nmina,nmaxa,typea,kwpa,vaa,detecta,pfa); 
%PickThreshold: Finds the left edge of the right-most mode
%
%   [th,thmax] = PickThreshold(x,nmin,nmax,type,pf,kw,va)
%
%   x      Data in which mode is to be located.
%   nmin   Smallest possible number of points within the mode.
%          Default = max(x).
%   nmax   Largest  possible number of points within the mode.   
%          Default = min(x);
%   type   Type of histogram: 1=MSE , 2=Matched filter (Default)
%   kw     Kernel width
%   va     Valley depth
%   detect Detect one major type of spikes(1) or all spikes(2)
%   pf     Plot flag argument: Default=0 (No figures)
%
%   th     Threshold selected by the algorithm
%   thmax  Maximum threshold limit
%
%   Is sensitive to outliers - picks critical paramters based on 
%   range of data 
%
%   Version 1.00 SK
%
%   See also DetectMaxima, DetectMinima, and PowerPeaks.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help PickThreshold;
    return;
    end;
   
%====================================================================
% Process Function Arguments
%====================================================================   
nmin = min(x);
if exist('nmina') & ~isempty(nmina),
    nmin = nmina;
    end;   
    
nmax = max(x);
if exist('nmaxa') & ~isempty(nmaxa),
    nmax = nmaxa;
    end;   

detect = 1;
if exist('detecta') & ~isempty(detecta),
    detect = detecta;
    end;   
    
type=2;                                                    % Default measure type = Crosscorrelation
if exist('typea') & ~isempty(typea),
    type = typea;
    end
            
kwp = 0.70;                                                % Kernel width parameter (fraction of IQR)
if exist('kwpa') & ~isempty(kwpa),
    kwp = kwpa;
    end

va = 0.5;                                                  % Depth of a valley between two modes
if exist('vaa') & ~isempty(vaa),
    va = vaa;
    end

pf = 0;                                                    % Default - no figures
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Designer-Specified Parameters
%====================================================================
nb  = 1000;                                                % No. of bins

%====================================================================
% Preprocessing
%====================================================================
if type==1,
    x = -x;
    end;

x   = x(:);                                                % Make x a column vector, if not already
nx  = length(x);                                           % Number of data points
xrg = max(x)-min(x);                                       % Range of the data
bw  = xrg/nb;                                              % Bin width
kw  = kwp*iqr(x);                                          % Actual kernel width

%====================================================================
% Setting threshold limits
%====================================================================
xs    = sort(-x);                                          % Sort in descending order
xs    = -xs;                                               % Convert back to original coordinates in ascending order  
thmax = xs(nmin);                                          % Maximum threshold
thmin = xs(min(nmax,length(xs)));                          % Minimum threshold
		
%====================================================================
% Main Decision Logic for Selecting a Threshold
%====================================================================
xt = min(x)-kw:bw:max(x)+kw;                               % Pick the bin centers
[f,xt] = KernelDensity(x,xt,kw);                           % Estimate the PDF with a kernel smoother
f   = f(:);                                                % Convert to column vector, if not already

pim   = DetectMaxima(f,[],0);                              % Locate local maxima
[jnk, spi] = sort(f(pim),1,'descend');
spimu = pim(spi);                                          % Sort maxima indices based on f values, unlimited
pip   = spimu(1);                                          % Primary mode peak index
ind   = find( xt(spimu) >= xt(pip) & xt(spimu) <= thmax);  % Figure out which spimu are within the valid threshold range
spim  = spimu(ind);                                        % Sorted peak indices to the right of the primary mode and less than the maximum threshold (limited)

pin = DetectMinima(f,[],0);                                % Locate local minima
pin = pin(xt(pin)>=thmin & xt(pin)<=thmax);
[jnk, spin] = sort(f(pin),1,'descend');
spin= pin(spin);                                           % Sort local minima based on f values
lpin = [1;pin;length(xt)];                                 % Include two edges of f into local minima

%=====================================================
% Select a threshold 
%=====================================================
if length(pim)==1,                                         % If there is only one maximum
    th = thmax;                                            % threshold equals threshold maximum
elseif length(pim)==2,                                     % If there are only two maxima
    thid = find(f(pim(1):pim(2))==min(f(pim(1):pim(2))));  % Find the minimum value between the maxima
    thid = pim(1)+thid-1;                                  % Shift index since was only a search in the range pim(1):pim(2)
    if thid>length(xt)
        thid=length(xt);
    end;
    if xt(thid)<=thmin,
        th = thmin;
    elseif xt(thid)>=thmax,
        th = thmax;
    else
        th = xt(thid);                                     % Set a threshold at the lowest point between them
    end;
elseif length(pim)>=3,                                     % If there are more than two maxima  
    c2 = 1;                                                % Counter for threshold candidates
    P  = [];
    for c1=1:length(spin),
        cpin = find(lpin==spin(c1));                       % Choose a local minimum according to its f value
        d    = pim-lpin(cpin(1,1));                        % Locate two local maxima closest to the chosen local minimum
        nd   = d(d<0);
        pd   = d(d>0);
        lmax = pim(d==max(nd));
        rmax = pim(d==min(pd));
       
        if f(lpin(cpin)) < va * min(f(lmax),f(rmax)),      % If the chosen local minimum's f value is less than va x the smaller peak value among two peaks
            P(c2) = spin(c1);
            c2    = c2 + 1;            
        else                                               % Otherwise,
            id   = find(f(pim)==min(f(lmax),f(rmax)));     % Locate the index of the lower of the two maxima (left & right)
            pim  = [pim(1:id-1);pim(id+1:end)];            % Eliminate this maximum from the list of maxima
            lpin = [lpin(1:cpin-1);lpin(cpin+1:end)];      % Eliminate the minima from the list
        end;        
    end;

    if isempty(P),                                         % If P is empty, th = thmax;
        th = thmax;
    else    
        if detect == 1,
            th = xt(max(P));                                   % Detect the main spikes
        else
            th = xt(min(P));                                   % Detect all spikes
        end
    end;
end
 
%====================================================================
% Plot How the Algorithm Makes its Decision
%====================================================================    
if pf,
    figure;
    FigureSet;
    
    for c1=1:2,
        bcp = min(x):bw*10:max(x);                             % Pick the bin centers
        [bhp,bcp] = hist(x,bcp);                               % Generate histogram
        bhp = bhp/(sum(bhp)*bw*10);                            % Normalize histogram to have unit area
        switch c1,
        case 1,
            ax = axes('Position',[0.10 0.79 0.85 0.20]);
            ymax = max(bhp)*1.05;                                  % Y limit for upper axis
        case 2,
            ax = axes('Position',[0.10 0.10 0.85 0.59]);
            if length(spimu)>1,
                ymax = 1.50*max(f(spimu(2:end)));
            else
                ymax = max(bhp)*1.05; 
                end;
            end;
        
        h = bar(bcp,bhp,1);
        set(h,'FaceColor',0.7*[1 1 1]);                        % Light gray 
        set(h,'EdgeColor',0.6*[1 1 1]);                        % Light gray    
        lc = 0;

        hold on;        
            h = plot(thmin*[1 1],[0 ymax],'k-.',thmax*[1 1],[0 ymax],'k-.'); % Threshold Limits
            set(h,'LineWidth',1.5);
            ha = [h(1)];
            lc = lc +1;
            lgl{lc} = 'Threshold Boundaries';                  % Legend Label

            h = plot(xt,f,'k');                                % Estimated PDF
            set(h,'LineWidth',2.5);
            set(h,'Color', 0.1*[1 1 1]);
            ha = [ha;h];
            lc = lc +1;
            lgl{lc} = 'Estimated PDF';                         % Legend Label

            h = plot(xt(spim),f(spim),'ks');                   % All maxima
            set(h,'MarkerFaceColor',[1.0 0.9 0.9]);
            set(h,'MarkerEdgeColor',[0.5 0.0 0.0]);
            set(h,'MarkerSize',6);
            ha = [ha;h];
            lc = lc +1;
            lgl{lc} = 'Maxima';
            
            [jnk,id] = min(abs(th-xt));                 
            h = plot([th th],[0 ymax],'g');
            set(h,'LineWidth',1.2);
            set(h,'Color',[0.0 0.5 0.0]);
            ha = [ha;h(1)];      
            lc = lc + 1;
            lgl{lc} = 'Threshold';
            
            if ~isempty(pin),
                h = plot(xt(pin),f(pin),'ko');                 % Local minima = Threshold candidates
                set(h,'MarkerFaceColor',[0.9 1.0 0.9]);
                set(h,'MarkerEdgeColor',[0.0 0.5 0.0]);
                set(h,'MarkerSize',6);
                ha = [ha;h];
                lc = lc +1;
                lgl{lc} = 'Threshold candidates';
                end

        hold off;        

        xlim([min(bcp) max(bcp)]);
        ylim([0 ymax]);
        %FormatTicks('%1.0f','%1.2f');
        switch c1,
        case 1,
        case 2,
            legend(ha,lgl,'Location','NorthWest');            
            xlabel('Cross-correlation (Scaled)');
            ylabel('Freqeuncy (Normalized)');
            end;
        %set(gca,'layer','top');
        AxisSet(8);    
        box off;  
        end;
    
%     if type==2,                                            % Plot a small window to show the detail
%         h=axes;
%         set(h,'position',[0.55 0.3 0.2 0.2]);
%         h = bar(bcp(ceil(length(bcp)*0.5):ceil(length(bcp)*0.85)),bhp(ceil(length(bhp)*0.5):ceil(length(bhp)*0.85)),1);
%         set(h,'FaceColor',0.7*[1 1 1]);                    % Light gray 
%         set(h,'EdgeColor',0.6*[1 1 1]);                    % Light gray 
% 
%         hold on;
%             h = plot(xt(ceil(length(xt)*0.5):ceil(length(xt)*0.85)),f(ceil(length(f)*0.5):ceil(length(f)*0.85)),'k'); % Estimated PDF    
%             set(h,'LineWidth',2);
%             set(h,'Color', 0.1*[1 1 1]);
%             g = plot(xt(pin),f(pin),'ko');                 % Local minima = Threshold candidates
%             set(g,'MarkerFaceColor',[0.9 1.0 0.9]);
%             set(g,'MarkerEdgeColor',[0.0 0.5 0.0]);
%             set(g,'MarkerSize',6);
% %             [h] = VArrow(f(id)+max(f)*.2, f(id)*1.2, th,[],max(xt)*0.01,ymax*0.05); % Threshold arrow
% %             [h] = VArrow(f(id)*50, f(id)*1.2, th, [],max(xt)*0.02,ymax*0.015); % Threshold arrow
%         hold off;
%         xlim([xt(ceil(length(xt)*0.5)) xt(ceil(length(xt)*0.85))]);
%         ylim([0 max(f(ceil(length(f)*0.5):ceil(length(f)*0.85)))]);
%         set(gca,'layer','top');
%         AxisSet(8);    
%         zoom on;
%         box off;    
%     end    
    end
  
%====================================================================
% Postprocessing
%====================================================================
if type==1,
    th = -th;
    end;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('th');
end;