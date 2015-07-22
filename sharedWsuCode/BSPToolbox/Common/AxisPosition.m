function [ha] = AxisPosition(opa)
%AxisPosition: Creates or positions axes in common configurations
%
%   [ha] = AxisPosition(op)
%
%   op   Type of arrangement: 1=Single axis (default).
%
%   ha   Array of axes that were positioned within the figure,
%        in order of left to right then top to bottom.
%
%   Arranges different MATLAB axes within a figure to more fully
%   occupy the figure area. This reduces much of the white space
%   that results when using MATLAB's default values.
%
%   1   The default position for a single axis uses only 63% of
%       the figure area. This positions a single axis to use
%       is 78.3% while still leaving room for the axis labels.
%
%   Example: Generate a generic plot for use in a PowerPoint slide 
%   and apply AxisPosition.
%
%      load ICP.mat; 
%      k = 1:round(fs)*5;
%      FigureSet(1,'PPT'); 
%      plot((k-1)/fs,icp(k));
%      xlabel('Time (sec)');
%      ylabel('ICP (mmHg)');
%      title('ICP Signal Segment (5 sec)');
%      AxisPosition;
%      AxisSet;
%      print -dtiff -r150 Generic; 
%
%   Edward Tufte, "The Visual Display of Quantitative Information, 
%   Cheshire", CT: Graphics Press, 1997.
%
%   Version 0.00.00.23 JM
%
%   See also SUBPLOT, FigureSet and AxisSet.

%====================================================================
% Process function arguments
%====================================================================
op = 1;
if exist('opa') & ~isempty(opa),
	op = opa;
	end;

%====================================================================
% Main Routine
%====================================================================
switch op
case 1,
    ha = gca;
    if length(ha)~=1,
        clf;
        ha = axis;
        end;
    set(ha,'Position',[0.08 0.08 0.90 0.87]);
otherwise
    error('First input argument was invalid.');
    end;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('ha');
    end;
