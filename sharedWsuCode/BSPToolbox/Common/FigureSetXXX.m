function [] = FigureSet(qda,wda,hta)
%FigureSet: Scales and positions figures in common configurations
%
%   [] = FigureSet(qd,wd,ht)
%
%   qd   Screen quadrant: 1=Top left (default), 2=Top Right, 
%        3=Bottom left, 4=Bottom right
%   wd   Width of figure (inches). Optionally 'PPT', 'wide','PDF',
%        or 'Slides'. Default=6.5 in.
%   ht   Height of figure (inches). Default=4.017 in.
%
%   Moves figure to specified quadrant. If only height or width is 
%   specified, the figure is scaled so that the ratio of the height 
%   to the width is the golden ratio. This is supposed to be the most 
%   visually pleasing. The optional automatic scaling for the wd 
%   argument are as follows. All of these set the plot orientation to 
%   portrait.
%
%   'Wide'    The figure occupies the entire top half of the screen 
%             and the ht to wd ratio is set to half the golden ratio.
%   'PPT'     The width is set to 8.5 and the height is set to 4.9. 
%             This works well for PowerPoint slides printed in tiff 
%             format.
%   'PDF'     This works well for printing a plot across an entire 
%             letter-sized page as postscript and then converting to 
%             PDF. Requires the width=11 and height=8.5 in Adobe 
%             Acrobat distiller (Settings>Job Options>Advanced).
%   'Slides'  Works well for LaTeX slides. Width = 4.5, Height = 2.7.
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
%   Version 1.00.01.10 JM
%
%   See also AxisSet and FormatTicks.

%====================================================================
% Preprocessing
%====================================================================
%orient portrait;
GoldenRatio = (sqrt(5)-1)/2; %0.6180

%====================================================================
% Process function arguments
%====================================================================
qd = 1;
if exist('qda') & ~isempty(qda),
	qd = qda;
	end;

wf  = 0;                                                   % Wide plot flag
wpf = 0;                                                   % Wide paper format flag
wd = 6.50;
if exist('wda') & isstr(wda),
    if strcmpi(wda,'PPT'),                                 % PowerPoint
        wd = 8.5;
        ht = 5.17;
    elseif strcmpi(wda,'Wide'),                            % Full width of screen
        wf = 1;
        ht = wd*GoldenRatio/2;
    elseif strcmpi(wda,'PDF'),                             % PDF
        wd  = 10.5;
        ht  = 8;
        wpf = 1;
    elseif strcmpi(wda,'LTX') | strcmpi(wda,'Slides'),     % LaTeX Seminar Slides
        wd = 4.5;
        ht = 2.6;        
        end;
else
    if exist('wda') & ~isempty(wda),
	    wd = wda;
	    end;
    ht = wd*GoldenRatio;
    if exist('hta') & ~isempty(hta),
	    ht = hta;
	    end;
    end;

%====================================================================
% Determine Placement of Figure on the Screen
%====================================================================    
ss  = get(0,'ScreenSize');
if wd>ht,
    if wf,
    	swd  = ss(3) - 10;     % Screen width
	    sht  = swd * ht/wd;    % Screen height 
    else
    	swd  = ss(3)/2 - 10;
	    sht  = swd * ht/wd;
        end;
else
	sht  = ss(4)/2;
	swd  = sht * wd/ht;
	end;

drawnow;    
fh = get(0,'Children');  % Figure handles
fh = fh(fh~=gcf);        % Exclude the current figure being worked on 
nf = length(fh);         % No. figures
fp = get(fh,'Position'); % Figure positions
switch qd
	case 1,
		le  = 10;
		be  = ss(4)-sht-70;
        c = 0;
        if nf>1,
            for c1 = 1:nf,   % Count the number of figures in this approximate position
                if fp{c1}(1)==le & fp{c1}(2)>ss(4)*0.25,
                    c = c + 1;
                    end;
                end;
        elseif nf==1 & fp(1)==le & fp<ss(4)*0.75,
            c = 1;
            end;
        be = be - c*10;
	case 2,
		le = ss(3)/2+8;
		be = ss(4)-sht-70;
        c = 0;
        if nf>1,
            for c1 = 1:nf,   % Count the number of figures in this approximate position
                if fp{c1}(1)==le & fp{c1}(2)>ss(4)*0.25,
                    c = c + 1;
                    end;
                end;
        elseif nf==1 & fp(1)==le & fp<ss(4)*0.75,
            c = 1;
            end;
        be = be - c*10;    
	case 3,
		le = 10;
		be = ss(4)-2*(sht+95);
	case 4,
		le = ss(3)/2+8;
		be = ss(4)-2*(sht+95);
    end;
sp  = [le be swd sht]; % Put in upper left hand corner of screen

%drawnow;
set(gcf,'Position',sp);

%====================================================================
% Place the figure in the center of the paper
%====================================================================  
set(gcf,'PaperUnits','Inches');
le = 8.5/2 - wd/2;       % Left Edge
be = 11/2  - ht/2;       % Bottom Edge
if wpf, % Wide paper format
    le = 11.0/2 - wd/2;  % Left Edge
    be =  8.5/2 - ht/2;  % Bottom Edge
    end;
set(gcf,'PaperPosition',[le be wd ht]);
drawnow;