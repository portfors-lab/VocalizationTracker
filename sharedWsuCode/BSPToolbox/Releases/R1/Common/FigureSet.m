function [] = FigureSet(qda,wda,hta)
%FigureSet: Scales and positions figures in common configurations
%
%   [] = FigureSet(qd,wd,ht)
%
%   qd   Screen quadrant: 1=Top left (default), 2=Top Right, 
%        3=Bottom left, 4=Bottom right
%   wd   Width of figure (inches). Optionally 'PPT', 'wide', or 
%        'PDF'. Default=6.5 in.
%   ht   Height of figure (inches). Default=4.017 in.
%
%   Moves figure to specified quadrant. If only height or width is 
%   specified, the figure is scaled so that the ratio of the height 
%   to the width is the golden ratio. This is supposed to be the most 
%   visually pleasing. The optional automatic scaling for the wd 
%   argument are as follows. All of these set the plot orientation to 
%   portrait.
%
%   'Wide' The figure occupies the entire top half of the screen and 
%          the ht to wd ratio is set to half the golden ratio.
%   'PPT'  The height is set to 8.5 and the height is set to 4.9. 
%          This works well for PowerPoint slides printed in tiff 
%          format.
%   'PDF'  This works well for printing a plot across an entire 
%          letter-sized page as postscript and then converting to PDF. 
%          Requires the width=11 and height=8.5 in Adobe Acrobat
%          distiller (Settings>Job Options>Advanced).
%
%   The default width argument works well for tiff files embedded in 
%   MS Word. 
%
%   Example: Generate a generic plot for use in a PowerPoint slide 
%   and apply FigureSet.
%
%      load ICP.mat; 
%      k = 1:round(fs)*5;
%      FigureSet(1,'PPT'); 
%      plot((k-1)/fs,icp(k));
%      print -dtiff -r150 Generic; 
%
%   Edward Tufte, "The Visual Display of Quantitative Information, 
%   Cheshire", CT: Graphics Press, 1997.
%
%   Version 1.00 JM
%
%   See also AxisSet.
