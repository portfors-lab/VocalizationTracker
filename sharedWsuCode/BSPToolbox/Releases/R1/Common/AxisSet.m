function [] = AxisSet(fsa,fna,aha)
%AxisSet: Converts axes to a specified font size and type
%
%   [] = AxisSet(fs,fn,ah)
%
%   fs   Font size (pts). Default=10 pts.
%   fn   Font name (text). Default='Time New Roman'.
%   ah   Axis handle. Default=all axes in current figure. 
%
%   Converts all text to the specified font and font size. Handy for
%   generating converting all of the text on a figure to a more
%   professional font like Times New Roman. Font size must be 1 or
%   greater. The possible fonts that can be used are listed in the
%   property editor for axes under the style tab and listed as 
%   "Font name". 
%
%   If the legend command is used with this command, it should be 
%   called after AxisSet so that the legend box size will be scaled
%   appropriately.
%
%   Example: Generate a standard plot and adjust the axes to use 
%   10 pt Comic Sans MS font.
%
%      load ICP.mat; 
%      FigureSet;
%      k = 1:round(fs)*5;
%      plot((k-1)/fs,icp(k));
%      xlabel('Time (s)');
%      ylabel('ICP (mmHg)');
%      text(2,1,'Text');
%      AxisSet(8,'Comic Sans MS');
%      legend('ICP (mmHg)');
%
%   Edward Tufte, "The Visual Display of Quantitative Information, 
%   Cheshire", CT: Graphics Press, 1997.
%
%   Version 1.00 JM
%
%   See also FigureSet.
