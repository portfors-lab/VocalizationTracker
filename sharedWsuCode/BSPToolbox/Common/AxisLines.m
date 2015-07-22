function [h] = AxisLines(sta)
%AxisLines: Draws a vertical and horizontal lines through the origin.
%
%   [h] = AxisLines(st)
%
%   st   Style of lines. Default = 'k:' (Black dotted line).
%
%   h    Handles to line objects.
%
%   Plots lines that span the limits of the plot and go through the
%   origin. Useful for quickly showing the reference to the origin.
%   This should be called after the limits of the x and y axes have
%   been set.
%
%   Example: Generate a scatter plot of random variable and add axis 
%   lines to it that use a red dashed line for the style and a 
%   linewidth of 2.
%
%      figure;
%      FigureSet; 
%      plot(randn(50,1),randn(50,1),'.');
%      h = AxisLines('r--');
%      set(h,'LineWidth',2.0);
%
%   Desmond J. Higham and Nicholas J. Higham, "MATLAB Guide," 2000.
%
%   Version 1.00 JM
%
%   See also Circle, xlim, ylim, axis, and grid.

%====================================================================
% Process function arguments
%====================================================================
st = 'k:';    
if exist('sta') & ~isempty(sta),
    st = sta;
    end;
    
%====================================================================
% Plot
%====================================================================
hf = ishold; % Hold flag
if ~hf,
    hold on;
    end;
    
ha = get(gca,'Children');
h  = [];
if max(ylim)>0 & min(ylim)<0,
    h1 = plot(xlim,0*ylim,st);
    set(h1,'LineWidth',0.5);
    set(h1,'Clipping','Off');
    h  = h1;
    end;
if max(xlim)>0 & min(xlim)<0,    
    h2  = plot(0*xlim,ylim,st);
    set(h2,'LineWidth',0.5);
    set(h2,'Clipping','Off');
    h = [h;h2];
    end;
if max(zlim)>0 & min(zlim)<0,    
    h3  = plot3(0*xlim,0*ylim,zlim,st);
    set(h3,'LineWidth',0.5);
    set(h3,'Clipping','Off');
    h = [h;h3];
    end;    
    
set(gca,'Children',[ha;h]); % Put axis lines on bottom

if ~hf, 
	hold off;
	end;
    
%====================================================================
% Process Output Arguments
%====================================================================      
if nargout==0,
    clear h;
    return;
    end;    
