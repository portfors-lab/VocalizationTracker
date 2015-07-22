function [] = Chart(fn)
%Chart: Generates the pseudo-Gantt chart for project tracking
%
%   [] = Chart(fn);
%
%   fn   File name.
%
%   Generates a chart similar to a Gantt chart used for project 
%   tracking. The function reads the Excel spreadsheet specified
%   by the filename. This should be in the same format as the 
%   example file 'ManuscriptSchedule.xls' in the Documentation
%   subdirectory of the BSP toolbox.
%
%   Example: Generate a chart using the example schedule file.
%
%      Chart('ManuscriptSchedule.xls');
%
%   R. L. Harris, "Information Graphics. A comprehensive illustrated
%   reference," Oxford University Press, 1999.
%
%   Version 1.00 JM
%
%   See also xlsread and FigureSet.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help Chart;
    return;
    end;

%====================================================================
% Load the File
%====================================================================
[N,T] = xlsread(fn);                       % Load numeric (N) and text (T) data from file

al = T(1,1);                                               % Author list
mt = T(1,2);                                               % Manuscript title


T = T(3:end,:);                                            % Strip off the headers

if size(T,1)~=size(N,1) | size(T,2)~=11 | size(N,2)~=8,
    error('Spreadsheet formatted incorrectly.');
    end;    

%====================================================================
% Preprocessing
%====================================================================
ns = size(N,1);                                            % Number of steps
emax = max([N(:,3)-N(1,2);N(:,5)-N(1,2)]);                 % Maximum number of days used for project (xmax)

%====================================================================
% Draw the Chart
%====================================================================
figure;
FigureSet(1);
dt  = zeros(ns,1);
da  = zeros(ns,1);
sd  = N(1,2) ;                                           % Start date   (Excel time)
td  = now - (sd+datenum('30-Dec-1899'));                 % Today's date (days since start of project)
h   = plot([td td],[0 ns+1],'b');
set(h,'LineWidth',1.5);
ylim([0 ns+1]);
xlim([0 emax]);
hold on;

[y,m,d] = datevec(sd+datenum('30-Dec-1899'));
d = 1;                                                     % Advance to the first day
m = m +1;                                                  % of the next month

md = datenum(y,m,d)- datenum('30-Dec-1899') -sd;           % Beginning of the month date
while 1,
    h = plot([md md],ylim,'k:');
    set(h,'LineWidth',0.2);
    st = datestr(datenum(y,m,d),'mmmyy');
    h = text(md,0,st);
    set(h,'Rotation',90);
    set(h,'VerticalAlignment','Top');
    set(h,'HorizontalAlignment','Right');
    m = m + 1;
    md = datenum(y,m,d)- datenum('30-Dec-1899') -sd;
    if md>emax,
        break;
        end;
    end;
hold off;

for c1 = 1:ns,
    h = text(-0.025*emax,c1,T(c1,1));                      % Add text to vertical axis with brief description of this step
    set(h,'HorizontalAlignment','Right');              
    
    bg     = N(c1,2)-sd;                                   % Scheduled start of this task (days relative to start date)
    ed     = N(c1,3)-sd;                                   % Scheduled end   of this task (days relative to start date)
    dt(c1) = ed-bg;                                        % Scheduled task duration
    if ed>bg,                                              % If task has finite duration
        h = patch([bg ed ed bg],c1+[-0.45 -0.45 0.45 0.45],0.5*[1 1 1]);
    else                                                   % Otherwise task is a milestone - plot as a diamond
        mp = (bg+ed)/2;
        h = patch([mp mp+2 mp mp-2],c1+[-0.45 0 0.45 0],0*[1 1 1]);
        end;
    
    if ~isnan(N(c1,5)),                                    % If task has actually started,
        bg = N(c1,5) - sd;                                 % Actual start of this step (days relative to start date)
        if isnan(N(c1,6)),
            ed = now - datenum('30-Dec-1899') - sd;        % Pick now as the end date if task is not yet completed
        else
            ed = N(c1,6) - sd;                             % Actual end date of this step
            end;
        da(c1) = ed-bg;                                    % Actual task duration to date
        if da(c1)<=dt(c1),
            clr = [0 1 0];
        elseif da(c1)<=dt(c1)+7,
            clr = [1 1 0];
        else
            clr = [1 0 0];
            end;           
        if ed>bg,
            h = patch([bg ed ed bg],c1+[-0.25 -0.25 0.25 0.25],clr);
        else
            mp = (bg+ed)/2;
            h = patch([mp mp+2 mp mp-2],c1+[-0.25 0 0.25 0],clr);
            end;
        end;
    end;
box off;
xlabel('Days');
title(sprintf('"%s", %s',char(mt),char(al)));
set(gca,'YTick',1:ns);
set(gca,'Position',[0.15 0.05 0.98-0.15 0.98-0.05]);
set(gca,'YDir','Reverse');
orient landscape;
AxisSet;
    