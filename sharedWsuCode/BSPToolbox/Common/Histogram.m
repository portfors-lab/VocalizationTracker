function [p,u] = Histogram(x,bwa,kwa,pfa);
%Histogram: Generates an advanced histogram with estimated PDF.
%
%   [p,u] = Histogram(x,bw,kw,pf);
%
%   x     Input signal.
%   bw    Bin width of histogram. Default = IQR/20.
%   kw    Kernel smoother width. Default = bw*2.    
%   pf    Plot flag: 0=none(default), 1=new figure, 2=current figure.
% 
%   p     Kernel estimate of density.
%   u     Values at which estimate was made.
%
%   Returns the peaks in the local smoothed signal energy. Signal energy 
%   is defined here as the absolute value of the signal to the power of p. 
%   sorted vector x that are at least as big as xm (optional). Uses a 
%   Lowpass filter type 4 to prevent ringing after impulses.
%
%   Example: Find the energy peaks in an electrocardiogram signal.
%
%      load Tremor;
%      [pi,y] = PowerPeaks(x,fs,150);
%      Histogram(y(pi));
%
%   J. Simonoff, Smoothing Methods in Statistics. New York, 
%   NY: Springer-Verlag, 1996.
%
%   Version 1.00.22 JM
%
%   See also HIST, BAR, and SmoothSeries.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help Histogram;
    return;
    end;

nx = length(x);
if nx<5,
    error('Data vector contains fewer than 5 points.\n');
    return;
    end;

%====================================================================
% Process function arguments
%====================================================================
if exist('bwa') & ~isempty(bwa),
    bw = bwa;
else
    bw = diff(prctile(x,[25 75]))/20;                      % Default
    if bw==0,
        bw = (max(x)-min(x))/20;                           % In case the IQR==0
        end;
    if bw==0,                                              % In case max(x)==min(x)
        bw = 1;
        end;
    end;
    
kw = bw*2;   
if exist('kwa') & ~isempty(kwa),
    kw = kwa;
    end;
       
pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
  
%====================================================================
% Preprocessing
%==================================================================== 
iqr = diff(prctile(x,[25 75]));                            % Interquartile range
md  = median(x);                                           % Median
ul  = max(md - 1.5*iqr,min(x));                            % Lower boundary
uu  = min(md + 1.5*iqr,max(x));                            % Upper boundary
if ul==uu,
    ul = ul - bw;
    uu = uu + bw;
    end;
    
st  = bw/10;                                               % Step size
    
u  = ul:st:uu+bw/10;                                       % Points to estimate density at
uh = ul:bw:uu;                                             % Histogram bins
uh = [ul-bw,uh,uu+bw];
    
    
%====================================================================
% Estimate Density
%====================================================================
yh = hist(x,u);                                            % Put the data in sparse, sorted bins
yh = yh/(sum(yh)*st);                                      % Normalize to have unit area
p  = SmoothSeries(u,yh,u,kw);                              % Smooth the histogram
    
%====================================================================
% Plot Histogram
%====================================================================
if pf>0,
    if pf==1,
        figure;
        end;
    FigureSet;
        
    %----------------------------------------------------------------
    % Top Plot (wider range)
    %----------------------------------------------------------------    
    h = axes('Position',[0.10 0.91 0.88 0.05]);
    tl = ul - (uu-ul);
    tu = uu + (uu-ul);
    ut = tl:bw:tu;
    yu = hist(x,ut);
    yu = yu/(sum(yu)*bw); % Normalize to have unit area
    h  = bar(ut,yu,1);
    set(h,'FaceColor',0.0*[1 1 1]); % Black
    set(h,'EdgeColor',0.0*[1 1 1]); % Black    
    xlim([tl tu]);
    ylim([0 max(yu)*1.02]);
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    hold on;
        h = plot([ul ul],[0 max(ylim)],'r','LineWidth',1.2);
        h = plot([uu uu],[0 max(ylim)],'r','LineWidth',1.2);
        hold off;
    box off;
    
    %----------------------------------------------------------------
    % Labels on Top of the Top Plot
    %----------------------------------------------------------------        
    h = axes('Position',[0.10 0.95 0.88 0.05]);
    xlim([tl tu]);
    ylim([-1 1]);
    set(h,'Visible','Off');
    h = text(tl       ,0,sprintf('%4.2f',100*sum(       x<tl)/nx),'HorizontalALignment','Left');
    h = text((tl+ul)/2,0,sprintf('%4.2f',100*sum(tl<=x & x<ul)/nx),'HorizontalALignment','Center');
    h = text((ul+uu)/2,0,sprintf('%4.2f',100*sum(ul<=x & x<uu)/nx),'HorizontalALignment','Center');
    h = text((uu+tu)/2,0,sprintf('%4.2f',100*sum(uu<=x & x<tu)/nx),'HorizontalALignment','Center');
    h = text(tu       ,0,sprintf('%4.2f',100*sum(tu<=x       )/nx),'HorizontalALignment','Right');
    
    %----------------------------------------------------------------
    % Bottom Plot (main figure)
    %----------------------------------------------------------------
    h = axes('Position',[0.10 0.10 0.88 0.80]);
    yh = hist(x,uh);
    uh = uh(2:length(uh)-1);
    yh = yh(2:length(yh)-1);
    yh = yh/(sum(yh)*bw);                                  % Normalize to have unit area
    h  = bar(uh,yh,1);
    set(h,'FaceColor',0.7*[1 1 1]); % Light gray
    set(h,'EdgeColor',0.2*[1 1 1]); % Dark gray
    hold on;
        h = plot(u,p,'b');
        set(h,'LineWidth',1.5);
        hold off;
    box off;
    xlim([ul uu]);
    ylim([0 max(yh)*1.02]);
    ylabel('PDF Estimate');
    
    AxisSet;
    end    

%====================================================================
% Process Output Arguments
%====================================================================      
if nargout==0,
    clear('p','u');
    return;
    end;    
