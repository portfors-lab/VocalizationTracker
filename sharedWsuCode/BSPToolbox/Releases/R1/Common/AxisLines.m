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
