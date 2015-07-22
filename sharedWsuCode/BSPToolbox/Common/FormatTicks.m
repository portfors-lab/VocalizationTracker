function [] = FormatTicks(xfa,yfa,aha);
%FormatTicks: Converts ticks to specified format
%
%   [] = FormatTicks(xf,yf,ah);
%
%   xf   C-style format string for x-axis.
%   yf   C-style format string for x-axis.
%   ah   Axis handle. Default=all axes in current figure. 
%
%   Formats the x-axis and y-axis tick labels with the user-specified
%   format. This is a handy way of ensuring the same precision is 
%   used for each tick, which MATLAB does not do by default. 
%
%   This should be called only after the x- and y-axis limits have 
%   been set. If the axis is rescaled or the range of one of the axes
%   is changed after this function has been called, the tick labels 
%   may be incorrect. 
%
%   Example: Generate a standard plot and adjust the axes to use 
%   10 pt Comic Sans MS font.
%
%      load ICP.mat; 
%      FigureSet;
%      k = 5*fs:10*fs;
%      plot((k-1)/fs,icp(k));
%      xlabel('Time (s)');
%      ylabel('ICP (mmHg)');
%      xlim([5 10]);
%      ylim([9.5 13]);
%      box off;
%      FormatTicks('%4.1f','%4.1f');
%      AxisSet;
%
%   B. W. Kernighan and D. M. Ritchie, "The C programming language," 
%   Prentice Hall, 1988.
%
%   Version 1.00 JM
%
%   See also sprintf, AxisSet, and FigureSet.

%====================================================================
% Process function arguments
%====================================================================
if nargin<1,
    help FormatTicks;
    return;
    end;

xf = ''; % Default x-axis tick format = no change
if exist('xfa') & ~isempty(xfa),
	xf = xfa;
	end;
  
yf = ''; % Default y-axis tick format = no change
if exist('yfa') & ~isempty(yfa),
	yf = yfa;
	end;
	
ah = get(gcf,'Children');
if exist('aha') & ~isempty(aha),
    ah = aha;
    end;
    
%====================================================================
% Main Loop
%====================================================================   
for c1 = 1:length(ah),
    if ~isempty(xf),
        xt  = get(ah(c1),'XTick');
        xtl = {};
        for c2 = 1:length(xt),
            xtl(c2) = cellstr(sprintf(xf,xt(c2)));
            end;
        set(ah(c1),'XTickLabel',xtl);
        set(ah(c1),'XTickMode','Manual'); % Prevents new tick generation when printing
        end;
    if ~isempty(yf),
        yt  = get(ah(c1),'YTick');
        ytl = {};
        for c2 = 1:length(yt),
            ytl(c2) = cellstr(sprintf(yf,yt(c2)));
            end;
        set(ah(c1),'YTickLabel',ytl);
        set(ah(c1),'YTickMode','Manual'); % Prevents new tick generation when printing        
        end;
    end;