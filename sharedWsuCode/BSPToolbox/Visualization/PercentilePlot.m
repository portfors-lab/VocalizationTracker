function h = PercentilePlot(x,Y,pca,wf);
%PercenitlePlot: A percentile-based replacement for the box plot
%
%   [h] = PercentilePlot(x,Y,pc,wf);
%
%   x    Abscissa.
%   Y    Data with one column for each x.
%   pc   Percentiles to plot
%   wf   Widest bar width. Expressed as a fraction of the 
%        smallest inter-abscissa spacing that the wides. 
%        Default = 40. Should not exceed 50.
%
%   h    Vector of handles to plotted objects.
%
%   Generates a plot similar to a box plot, but based on percentiles.
%
%   Example: Plot the percentile plots of tremor segments.
%
%      ns = 20;
%      load Tremor.mat;
%      nx = length(x);
%      sl = floor(nx/ns);
%      Y = zeros(sl,ns);
%      t = zeros(ns,1);
%      i1 = 0;
%      for c1=1:ns, 
%          i0 = i1+1; 
%          i1 = i0+(sl-1); 
%          Y(:,c1) = x(i0:i1); 
%          t(c1) = ((i0+i1)/2-1)/fs;
%          end;
%      figure;
%      FigureSet;
%      PercentilePlot(t,Y,[50 75 90 95 99]);
%
%   R.L. Harris, "Information Graphics: A Comprehensive Illustrated 
%   Reference," Oxford University Press, 1999.
%
%   Version 1.00 JM
%
%   See also boxplot, AxisSet, AxisLines, and FigureSet.

% clear all;
% close all;
% FigureSet(1);       
% ns = 20;
% load Tremor.mat;
% nx = length(x);
% sl = floor(nx/ns);
% Y = zeros(sl,ns);
% t = zeros(ns,1);
% i1 = 0;
% for c1=1:10, 
%      i0 = i1+1; 
%      i1 = i0+(sl-1); 
%      Y(:,c1) = x(i0:i1); 
%      t(c1) = ((i0+i1)/2-1)/fs;
%      end;
% x = t;     
% pca = [50 75 90 95 99];

%====================================================================
% Error Checking
%====================================================================    
if nargin<2,
    help PercentilePlot;
    return;
    end;
    
if length(x)==0,
    error('Abscissa is empty.');
    end;

%====================================================================
% Process Function Arguments
%====================================================================
pc = [50 75 95 99];
if exist('pca') & ~isempty(pca),
    pc = pca;
    end;    
pc = sort(unique(pc));
    
bw = 40;
if exist('bwa') & ~isempty(bwa),
    bw = bwa;
    end;    
bw = min(bw,50);
bw = max(bw,1);

%====================================================================
% Preprocessing
%====================================================================
xs = min(diff(sort(unique(x))));                           % Minimum spacing between abscissa
ps = sort(pc);
nx = length(x);

%====================================================================
% Main Loop
%====================================================================
h = zeros(nx,1);
for c1=1:nx,
    for c2=1:length(pc),
        lp(c2) = prctile(Y(:,c1),50-pc(c2)/2);             % Lower ordinate for this percentile
        up(c2) = prctile(Y(:,c1),50+pc(c2)/2);             % Upper ordinate for this percentile
        if c2==1,                                          % If base percentage,
            wd  = xs*bw/100;                               % Set the width to 25% of the smallest abscissa spacing
            app = wd*(up(c2)-lp(c2))/pc(c2);               % Area per percentage
            if up(c2)==lp(c2),
                wd = 0;
                end;
        else
            pa = pc(c2)-pc(c2-1);                          % Percentage area
            if up(c2)==up(c2-1) && lp(c2)==lp(c2-1),
                wd = 0;
            else
                wd = pa*app/(up(c2)-up(c2-1)+lp(c2-1)-lp(c2)); % Width of bar  
                end;
            end;
        xw(c2) = wd;                                       % X-widths
        end; 
    xwr  = xw(end:-1:1);                                   % X-abscissa reversed
    lpr  = lp(end:-1:1);                                   % Lower percentiles reversed  
    upr  = up(end:-1:1);                                   % Upper percentiles reversed
    id1  = floor(1:0.5:length(pc));
    id2  = floor(1.5:0.5:length(pc)+0.5);    
    hp = patch(x(c1)+[-xwr(id1) -xw(id2) xwr(id1) xw(id2)],[lpr(id2) up(id1) upr(id2) lp(id1)],'r');
    h(c1) = hp;
    %fprintf('Pausing...\n'); pause;
    end;
set(h,'FaceColor',0.3*[1 1 1]);
set(h,'EdgeColor',0.3*[1 1 1]);
set(h,'LineWidth',0.001);


    