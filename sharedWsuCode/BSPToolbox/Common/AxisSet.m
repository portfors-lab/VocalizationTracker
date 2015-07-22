function [] = AxisSet(fsa,fna,aha)
%AxisSet: Converts axes to a specified font size and type
%
%   [] = AxisSet(fs,fn,ah)
%
%   fs   Font size (pts). Default=10 pts.
%   fn   Font name (text). Default='Times New Roman'.
%   ah   Axis handle. Default=all axes in current figure. 
%
%   Converts all text to the specified font and font size. Handy for
%   converting all of the text on a figure to a more professional font 
%   like Times New Roman. Font size must be 1 or greater. The possible 
%   fonts that can be used are listed in the property editor for axes 
%   under the style tab and listed as "Font name". 
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
%   Version 1.01 JM
%
%   See also FormatTicks, FigureSet, and figureLatex.

%   Changelog
%   Date    Author  Description
%   ------------------------------------------------------
%   3-13-06 SP      No longer changes interpreter to LaTeX
%                   Fixed typo in help text

%====================================================================
% Process function arguments
%====================================================================
fs = 10; % Default font size
if exist('fsa') & ~isempty(fsa),
	if fsa<1,
		help AxisSet; 
        return;
	else
  	    fs = fsa;
		end;
	end;

fn = 'Times New Roman'; % Default font
if exist('fna') & ~isempty(fna),
	fn = fna;
	end;

ax = get(gcf,'Children'); % Default vector of handles to axes
if exist('aha') & ~isempty(aha), 
    ax = aha;
    end;
    
%====================================================================
% Main Loop
%====================================================================    
for cnt = 1:length(ax),
    ha = ax(cnt); % Current axis handle
    
    %disp(get(ax(cnt),'Type'));
    if ~strcmpi(get(ax(cnt),'Type'),'axes'),
        continue;
        end;
    
	set(ha,'FontUnits','points');
	set(ha,'FontSize',fs); 
	set(ha,'FontName',fn);
    if ~strcmpi(get(ha,'FontName'),fn),
        error('Invalid font name.');
        end;
   
	h = get(ha,'Xlabel');
	set(h,'FontUnits','points');
	set(h,'FontSize',fs);
	set(h,'FontName',fn);
    set(h,'VerticalAlignment','middle');
%    set(h,'Interpreter','LaTeX');
   
	h = get(ha,'Ylabel');
	set(h,'FontUnits','points');   
	set(h,'FontSize',fs);
	set(h,'FontName',fn);
%    set(h,'Interpreter','LaTeX');
    
	h = get(ha,'Zlabel');
	set(h,'FontUnits','points');   
	set(h,'FontSize',fs);
	set(h,'FontName',fn);
%    set(h,'Interpreter','LaTeX');
   
	h = get(ha,'Title' );
	set(h,'FontUnits','points');   
	set(h,'FontSize',fs);
	set(h,'FontName',fn);
    set(h,'VerticalAlignment','middle');
%    set(h,'Interpreter','LaTeX');
   
%	h = findobj(ha,'Interpreter','LaTeX');
    h = findobj(ha,'FontName',fn);
	set(h,'FontUnits','points');
	set(h,'FontSize',fs);
	set(h,'FontName',fn);           
	end;
