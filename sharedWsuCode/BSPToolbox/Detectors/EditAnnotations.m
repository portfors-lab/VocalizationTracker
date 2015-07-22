function [ei] = EditAnnotations(x,fs,dia,oia);
%EditAnnotations: Edit the annotations of detected features.
%   
%   [ei] = EditAnnotations(x,fs,di,oi);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   di   Detection indices (samples). Default = [].
%   oi   Other detection indices (cell array) that should be 
%        displayed.
%
%   ei   Edited indices.
%
%   A handly tool for manually editing the annotations of automatic
%   detectors. 
%
%   Example: Edit detected beats of a noisy electrocardiogram. 
%
%      load NoisyECG.mat; 
%      ei = EditAnnotations(ecg,fs,di);
%
%   Task Force of the European Society of Cardiology and the American 
%   Society of Pacing and Electrophysiology, "Heart rate variability:
%   Standards of measurement, physiological interpretation, and 
%   clinical use," Circulation, vol. 93, no. 5, pp. 1043-1065, Aug. 
%   1996.
%
%   Version  0.00.01.23 JM
%
%   See also ManualDetect, DetectionPlots, and Detectors. 

%=====================================================================
% Change Log
%=====================================================================
%
%   Author Date     :  Description
%   ------------------------------------------------------------------
%   DT 10.20.02: Added option t,T for zooming to an window whose
%                length is approximately that of an average 
%                interval, Added lines 323-28
%   DT 10.20.02: Added try and catch blocks around entire user 
%                interface code to output data in case of any
%                error. added lines 106-8, 338-53
%   DT 10.20.02: Added arrow options for symmetric zoom, unzoom,
%                right and left. change lines 299,303,313,318 
%   JM 10.29.02: Added more comments, changed l and r commands to
%                move a full window length, cleaned up t command,
%                added expected error boundaries for rate and
%                interval plots, eliminated bugs in bottom plotting
%                routines.

%=====================================================================
% Process function arguments
%=====================================================================
if nargin<1,
    help EditAnnotations;
    return;
    end;

ei = [];
if exist('dia') & ~isempty(dia),
    ei = dia;
    ei = ei(:); % Force to be a column vector
    end;

oi = {};
if exist('oia') & ~isempty(oia),
    oi = oia;
    end;
    
if nargout==0,
    warning('No output argument specified. Annotations will not be saved.');
    end;
    
%=====================================================================
% Author-specified Parameters
%=====================================================================
mp  = 100e3;                                               % Maximum number of points allowed on a plot (speeds up plotting)
bp  = 1;                                                   % Bottom plot
sp  = 0;                                                   % Second plot
ntp = 100;                                                 % Maximum number of points to use to build template

%=====================================================================
% Preprocessing
%=====================================================================
ei = sort(ei);                                             % Ensure the indices are in increasing zorder

fvr = 0;                                                   % Initial value for "freeze vertical range" flag
sph = 0;                                                   % Second plot handle (zero if doesn't exist)
x   = x(:);                                                % Make into a column vector
nx  = length(x);
dc  = DistinctColors;

k  = (2:nx-1)';
%pi = k((x(k)> x(k-1) & x(k)>=x(k+1)) | (x(k)>=x(k-1) & x(k)> x(k+1))); % Indices of peaks and plateaus 
pi = DetectMaxima(x,[],1);

ts = 0;                                                    % Start of plot (seconds)
te = (nx-1)/fs;                                            % End   of plot (seconds)
ks = floor(ts*fs + 1);                                     % Start of plot (samples)
ke = ceil (te*fs + 1);                                     % End   of plot (samples)

%=====================================================================
% Main Loop 
%=====================================================================
fh = figure;
FigureSet(1,'Wide');
ha1 = axes('Position',[0.04 0.19 0.94 0.79]);
ha2 = axes('Position',[0.04 0.08 0.94 0.10]);

fprintf('\n');
fprintf('Left  : Add annotation & snap to nearest peak\n');
fprintf('Middle: Add annotation\n'); 
fprintf('Right : Remove annotation\n');
fprintf('1     : Display inter-beat intervals in bottom plot\n');
fprintf('2     : Display instantaneous rate in bottom plot\n');
fprintf('3     : Display peak amplitudes in bottom plot\n');
fprintf('4     : Display detection plot in another window\n');
fprintf('0     : Remove second plot\n');
fprintf('I,^   : Symmetric zoom in by a factor of 2\n');  
fprintf('A,v   : Symmetric zoom out by a factor of 2\n');
fprintf('Z     : Zoom in towards cursor by a factor of 2\n');
fprintf('T     : Zoom to a mean interval length window\n');
fprintf('L,<   : Move window left\n');
fprintf('R,>   : Move window right\n');
fprintf('F     : Freeze the vertical range at its current level\n'); 
fprintf('s     : Save annotations to temp.mat\n');
fprintf('x     : Exit function\n');

df = 0; % Done flag

%try,
while ~df,
    %-----------------------------------------------------------------
    % Plot Segment
    %-----------------------------------------------------------------
    axes(ha1);
    ss  = ceil((ke-ks+1)/100e3);                           % Step size - don't plot more than 100k points
    k   = ks:ss:ke;                                        % Index of samples to plot  
    t   = (k-1)/fs;                                        % Times of points plotted
    id  = find(ei>=ks & ei<=ke);                           % Indices of annotions that are within the plot range
    if length(id)==0,                                      % If no annotations in plot range, pick one
        [jnk,id] = min(abs(ei-(ks+ke)/2));                 % closest to midpoint
        end;
    if ~isempty(id) & id(1)>2,                             % If possible, prepend two annotations to list 
        id = [id(1)-2;id(1)-1;id];                         % of annoatation indices
        end;
    ni = length(id);        
    if ~isempty(id) & id(ni)<length(ei)-1,                 % If possible, postpend two annotations
        id = [id;id(ni)+1;id(ni)+2];                       % to list of annotation indices
        end;
    spi = ei(id);                                          % Segment peak indices - selected annotations
    spt = (spi-1)/fs;                                      % Segment plot times 
        
    ni  = length(spi);                                     % Number of annotations in plot range
    yl  = min(x(k));                                       % Maximum signal value within plot range
    yu  = max(x(k));                                       % Minimum signal value within plot range
    if yl==yu,                                             % If signal is constant,
        yl = yl*0.95 - 10*realmin;                         % pick lower & upper limits that insure
        yu = yu*1.05 + 10*realmin;                         % the vertical plot range is >0
        end;        
    yr  = yu-yl;                                           % Range of signal 
    yl  = yl - 0.05*yr;                                    % Add a small buffer to the lower limit
    yu  = yu + 0.05*yr;                                    % Repeat for the upper limit 
    if ~fvr,                                               % If user did not ask to keep vertical range frozen
        yrg = [yl yu];                                     % Update the vertical range
        end;
    
    cla;
    if ni>2 & ni<100,                                      % If the No. annotations is >2 & <100, 
        h  = plot(ones(2,1)*spt.',yrg'*ones(1,ni),'k:',... % Plot vertical dashed lines at annotation points
                  t,x(k),'b');              
    else                                                   % Otherwise, 
        h  = plot(t,x(k),'b');
        end;
    hold on;    
    
    for c1=1:length(oi),
        ri   = oi{c1};                                     % Reviewer indices
        id   = find(ri>=ks & ri<=ke);                      % Indices of annotions that are within the plot range
        ospi = ri(id);                                     % Other segment peak indices
        ospt = (ospi-1)/fs;                                % Other segmented peak times
        h = plot(ospt,x(ospi),'ro'); 
        set(h,'MarkerSize',(length(oi)-c1)*2+7);
        set(h,'Color',dc(c1+2,:));
        set(h,'MarkerFaceColor',dc(c1+2,:));
        end;
    
    if ni>2 & ni<100,                                      % If the No. annotations is >2 & <100, 
        h  = plot(spt,x(spi),'r.');              
    else                                                   % Otherwise, 
        h  = plot(spt,x(spi),'r.');
        end;
    if length(k)<500 & ni>2 & ni<100,                      % If the no. of signal points < 500, plot dots
        set(h,'Marker','.');
        end;
    hold off;
        
    xlim([t(1) t(length(t))]);                             % X-axis limits
    ylim(yrg);                                             % Y-axis limits 
    ylabel('Signal');
    set(ha1,'XTickLabel',[]);                              % Eliminate x-axis tick labels - is on bottom plot
    set(ha1,'TickLength',[0.002 0.0250])                   % Reduce the tick length 
    box off;
    
    axes(ha2);
    switch bp,
    %................................................................
    case 1,                                                % Plot intervals in bottom plot
    %................................................................
        ni = length(spi);
        if ni>1,
            ii = (spi(1:ni-1) + spi(2:ni))/2;              % Interval indices
            it = (ii-1)/fs;
            ix = diff(spi)/fs;                             % Interval times 
            mi = median(ix);                               % Median interval           
            yl = 0.65*mi;
            yu = 1.60*mi;
            if ni<100,
                h = plot([t(1) t(length(t))],0.75*mi*[1 1],'r:',...
                         [t(1) t(length(t))],1.50*mi*[1 1],'r:',...
                         ones(2,1)*spt.',[yl;yu]*ones(1,ni),'k:',...
                         it,ix,'g');
                nh = length(h);
                set(h(nh),'Marker','.');  
            else        
                h = plot([t(1) t(length(t))],0.75*mi*[1 1],'r:',...
                         [t(1) t(length(t))],1.50*mi*[1 1],'r:',...
                         it,ix,'g'); 
                end;
            ylim([yl yu]);
        else
            cla;
            end;
        ylabel('Interval (s)');
    %................................................................
    case 2,                                                % Plot instantaneous rate in bottom plot
    %................................................................
        ni = length(spi);
        if ni>1,
            ii = (spi(1:ni-1) + spi(2:ni))/2;              % Interval indices
            it = (ii-1)/fs;                                % Interval times
            ix = fs./diff(spi);                            % Instantaneous rate 
            mr = median(ix);                               % Median rate           
            yl = 0.65*mr;
            yu = 1.60*mr;
            if ni<100,
                h = plot([t(1) t(length(t))],0.75*mr*[1 1],'r:',...
                         [t(1) t(length(t))],1.50*mr*[1 1],'r:',...
                         ones(2,1)*spt.',[yl;yu]*ones(1,ni),'k:',...
                         it,ix,'m');                
                nh = length(h);
                set(h(nh),'Marker','.');        
            else        
                h = plot([t(1) t(length(t))],0.75*mr*[1 1],'r:',...
                         [t(1) t(length(t))],1.50*mr*[1 1],'r:',...
                         it,ix,'m');  
                end;
            ylim([yl yu]);
        else
            cla;
            end;
        ylabel('Rate (Hz)');        
    %................................................................
    case 3, % Plot peak amplitudes in bottom plot
    %................................................................
        ni = length(spi);
        if ni>0,
            spa = x(spi);
            p  = prctile(spa,[25 85]);
            if p(2)==p(1),
                yl  = min(spa);
                yu  = max(spa);
                if yl==yu,
                    yl = yl*0.95 - 10*realmin;
                    yu = yu*1.05 + 10*realmin;
                    end;
                yr  = yu-yl;
                yl  = yl - 0.05*yr;
                yu  = yu + 0.05*yr;
            else
                yl = p(1) - 2.5*(p(2)-p(1));
                yu = p(2) + 2.5*(p(2)-p(1));    
                end;
            if ni<100,
                h  = plot(ones(2,1)*spt.',[yl;yu]*ones(1,ni),'k:',spt,spa,'r');
                nh = length(h);
                set(h(nh),'Marker','.');
            else        
                h   = plot(spt,spa,'r');
                end;    
            ylim([yl yu]);
        else
            cla;
            end;
        ylabel('Amplitude');
        end;
    xlabel('Time (sec)');
    xlim([t(1) t(length(t))]);    
    set(ha2,'TickLength',[0.002 0.0250])    

    %-----------------------------------------------------------------
    % Deal with the Second Plot, if Any
    %-----------------------------------------------------------------    
    if sp~=0,
        if sph==0,
            sph = figure;
            FigureSet(4);
            end;
        figure(sph);
        switch sp,
        %............................................................
        case 4,
        %............................................................
            DetectionPlots(x,fs,ei,[],0);
        %............................................................
        case 5,
        %............................................................
            if length(ei)>1,
                iei = ceil(mean(diff(sort(ei))));          % Mean inter-event interval
                eid = (-iei:iei)';
                t   = eid/fs;
                ne  = length(ei);
                I   = eid*ones(1,ne) + ones(length(eid),1)*ei';
                I(I<=0) = 1;
                I(I>nx) = nx;
                h   = plot(t,x(I),'k');             
                xlim([min(t) max(t)]);
                ylim([min(min(x(I))) max(max(x(I)))]);
                box off;
                xlabel('Time (s)');
                ylabel('Signal');
                AxisSet;
                end;         
        %............................................................
        case 6,
        %............................................................           
            if length(pi)~=0 & length(ei)~=0,
                xmin = min(x([pi;ei]));
                xmax = max(x([pi;ei]));
                bc   = linspace(xmin,xmax,50)';            % Bin centers
                bhp = hist(x(pi),bc)';                     % Bin heights for detected indices 
                bhe = hist(x(ei),bc)';                     % Bin heights for detected indices 
                bw = min(diff(sort(unique(bc))));          % Bin width
                h1 = bar(bc,bhp,1.0);
                set(h1,'FaceColor',0.7*[1 1 1]);
                set(h1,'EdgeColor',0.7*[1 1 1]);
                hold on;
                    h2 = bar(bc,bhe,1.0);
                    set(h2,'FaceColor',0.3*[1 0 0]);
                    set(h2,'EdgeColor',0.3*[1 0 0]);
                    hold off;
                box off;
                xlabel('Signal');
                ylabel('Frequency');
                xlim([xmin-bw xmax+bw]);
                ymax = min(2*max(bhe),1.02*max(bhp));
                ylim([0 ymax]);
                AxisSet;
                legend([h1,h2],'Peaks','Annotations','Location','NorthWest');
                end;                
            end;
    elseif sp==0 && sph~=0,
        close(sph);
        sph = 0;
        end;
    figure(fh);                                            % Make the main figure the current figure again
        
    %-----------------------------------------------------------------
    % Obtain User Input
    %-----------------------------------------------------------------
    [ut, mag, ui] = ginput(1);
    
    if isempty(ui),
        fprintf('Empty!\n');
        continue;
        end; 
    
    switch ui,
        case 1, % Left mouse button
            us = round(ut*fs + 1);          % Convert cursor units from sec to samples        
            [jnk,id] = min(abs(us-pi));     % Find the closest corresponding peak
            up = pi(id);                    % User-specified peak (index)
            if ~isempty(find(up==spi)),     % Ignore if is already a peak
                continue;
                end;
            [jnk,id] = min(abs(up-ei)); % Find the closest corresponding peak
            ne = length(ei);
            if ei(id)>up,               % Figure out where to place it to preserve ordering
                ei = [ei(1:id-1);up;ei(id:ne)];
            else
                ei = [ei(1:id);up;ei(id+1:ne)];
                end;
        case 2, % Middle mouse button
            us = round(ut*fs + 1);          % Convert cursor units from sec to samples        
            if ~isempty(find(us==spi)),     % Ignore if is already a peak
                continue;
                end;
            [jnk,id] = min(abs(us-ei));     % Find the closest corresponding peak
            ne = length(ei);
            if ei(id)>us,                   % Figure out where to place it to preserve ordering
                ei = [ei(1:id-1);us;ei(id:ne)];
            else
                ei = [ei(1:id);us;ei(id+1:ne)];
                end;      
        case 3, % Right mouse button
            us = round(ut*fs + 1);          % Convert cursor units from sec to samples  
            if isempty(ei),
                continue;
                end;
            [jnk,id] = min(abs(us-ei));     % Find the closest corresponding peak
            ne = length(ei);
            ei = [ei(1:id-1);ei(id+1:ne)];    % Remove peak
        case {'1'}
            bp = 1;
        case {'2'}
            bp = 2;
        case {'3'}
            bp = 3;
        case {'4'}
            sp = 4;           
        case {'5'}
            sp = 5;               
        case {'0'}
            sp = 0;            
        case {'i','I',30}
            wl = round((ke-ks)/2);
            ks = ks + floor(wl/2);
            ke = ks + wl;
        case {'f','F'}
            fvr = ~fvr;                                    % Toggle the "freeze vertical range" flag
        case {'a','A',31}
            mp = round((ks+ke)/2);
            wl = (ke-ks);
            ks = max(1,mp-wl);
            ke = min(nx,mp+wl);    
        case {'z','Z'}
            st = max(10,round((ke-ks)/4));
            us = round(ut*fs + 1); % Convert cursor units from sec to samples
            ks = max( 1,us-st);
            ke = min(nx,us+st);      
        case {'r','R',29}
            wl = ke-ks;
            ke = min(nx,ke+wl);    
            ks = max( 1,ke-wl);
        case {'l','L',28}
            wl = ke-ks;
            ks = max( 1,ks-wl);
            ke = min(nx,ks+wl);
        case {'t','T'}
            id = find(ei>=ks-10*fs & ei<=ke+10*fs);
            mi = median(diff(ei(id))); % Median interval (samples)
            hi = round(mi/2);          % Half interval
            hi = max(5,hi);            % Limit to no less than 5 samples 
            us = round(ut*fs + 1);     % Convert cursor units from sec to samples
            ks = max( 1,us-hi);
            ke = min(nx,us+hi);  
        case {'s','S'}
            save('Temp.mat','ei');          
        case {'x','X'}
            break;
%             st = input('Are you sure? (y,n)','s');
%             if strcmpi(st,'y'),
%                 df = 1;
%                 break;
%                 end;
        end;         
    end;
         
% catch, % Catch block
%     fprintf(lasterr);
%     if nargout==1,
%         fprintf('Pausing... (hit ctrl-C to break, any other key to continue.\n');
%         pause;
%         return;
%     else
%         fprintf('\n\t');
%         saveData = input(['An error was encountered.''Save data [y,n]? '],'s');                
%         if saveData == 'y' | saveData == 'Y'
%             fprintf('\n\t');
%             filename = input('Filename : ','s');
%             save(filename,'ei');
%             end;
%         clear('ei');
%         end;
%     end; % End try-catch block   
    
%=====================================================================
% Process Outputs
%=====================================================================
if nargout==0,
    clear('ei');
    end;  
