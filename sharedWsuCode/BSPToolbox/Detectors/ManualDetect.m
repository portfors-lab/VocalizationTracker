function [mp,ibi] = ManualDetect(x);
%ManualDetect: Manual Peak Detector 
%   
%   [mp,ibi] = ManualDetect(x)
%
%   x     Input signal            
%
%   mp    Manual Detected Peaks, samples
%   ibi   Interbeat interval, samples
%
%   ManualDetect can be used to perform manual annotation of 
%   wave-component on signals such as ICP, ABP, ECG. The
%   function calculates the interbeat intervals, the first, 
%   and the second derivative, and plots them to aid
%   the manual annotation in regions with high artifact. 
%
%   Example: Perform manual detection on an the ICP signal. 
%
%      load ICP; 
%      [mp,ibi] = ManualDetect(icp);
%
%   Version  1.00 MA
%
%   See also EditAnnotations, PressureDetect, ECGDetectQRS and 
%   ECGDetectRInterbeat,

%=====================================================================
% Process function arguments
%=====================================================================
if nargin<1 | nargin>1,
    help ManualDetect;
    return;
    end;

%=====================================================================
% Process Inputs
%=====================================================================
R = 1;
x  = interp(x, R);                               
LD = length(x);
k  = 1:LD;

%=====================================================================
% Plotting the Input Signal
%=====================================================================
figure;
FigureSet(1,'Wide')
AxisSet(7,'Comic Sans MS');

subplot(4,1,1);
AxisSet(7,'Comic Sans MS');
title('Manual Detection');
ylabel('Magnitude');

subplot(4,1,2);
AxisSet(7,'Comic Sans MS');
title('Interbeat Intervals (IBI) in samples')
ylabel('IBI');

subplot(4,1,3);
AxisSet(7,'Comic Sans MS');
title('First Derivative');
ylabel('dx/dt');

subplot(4,1,4);
AxisSet(7,'Comic Sans MS');
title('Second Derivative');
ylabel('d(dx/dt)/dt');


subplot(4,1,1);
AxisSet(7,'Comic Sans MS');
plot(k, x, 'b');
axis tight;
title('Manual Detection');
ylabel('Magnitude');
xlabel('Samples');
hold on;
T = axis;

%%=====================================================================   
% Calculating and plotting the first derivative
%=====================================================================
FD = diff(x);
fdt= 1:length(FD);
subplot(4, 1, 3);
AxisSet(7,'Comic Sans MS');
title('First Derivative');
plot(fdt, FD, 'b');
ylabel('dx/dt');
hold on;

%=====================================================================
% Calculating and plotting the second derivative
%=====================================================================
SD = diff(x, 2);
sdt= 1:length(SD);
subplot(4, 1, 4);
AxisSet(7,'Comic Sans MS');
title('Second Derivative');
plot(sdt, SD, 'b');
ylabel('d(dx/dt)/dt');
hold on;


%=====================================================================
% User Menu
%=====================================================================
zoom on; 
fprintf('\n\n\t\t======================== MANUAL DETECTOR =====================\n');
fprintf('\t\t=     Pausing for zoom. Zoom in then hit enter to continue   =\n');
fprintf('\t\t==============================================================\n\n');
pause;
zoom off;
clc

fprintf('\n\n\t\t======================== MANUAL DETECTOR =====================\n');
fprintf('\t\t=     Left   Button: Label Peak                              =\n');
fprintf('\t\t=     Right  Buttom: Remove Peak                             =\n');
fprintf('\t\t=     1            : Replot inter-beat intervals (x2)        =\n');
fprintf('\t\t=     2            : Rezoom                                  =\n');
fprintf('\t\t=     0            : Save Detection and Quit                 =\n');
fprintf('\t\t==============================================================\n');



%=====================================================================
% Find all the peaks in the data set
%=====================================================================
PI = detectmaxima(x);
PI = PI(:);
MP = [];


%=====================================================================
% Manual Detection
%=====================================================================
stop = 0;
while ~stop
    [mi, mag, button] = ginput(1);
  
    %-----------------------------------------------------------------
    % No button or number detected
    %-----------------------------------------------------------------
    if isempty(button)
        fprintf('NO BUTTON OR NUMBER DETECTED \n');
        
    %-----------------------------------------------------------------
    % Left mouse buttom
    %-----------------------------------------------------------------
    elseif button == 1                               
        
        %.............................................................
        % Data and Peaks
        %.............................................................
        [mag, ind] = min(abs(mi-PI));              
        pi = PI(ind);
        MP = [MP; pi];
        subplot(4,1,1)
            hold on;
            h = plot(k(pi), x(pi), 'r.');
            AxisSet(7,'Comic Sans MS');
            title('Manual Detection');
            ylabel('Magnitude');
            set(h, 'MarkerSize', 15);                
            grid on;
       %.............................................................
       % Interbeat Intervals
       %.............................................................
       MP = sort(MP);
       P2P = diff(MP);
       dif = P2P;
       p2pi= (MP(1:length(MP)-1) + MP(2:length(MP)))/2;
       xl = xlim;
            subplot(4,1,2)
            hold on;
            h=plot(p2pi, dif);
            AxisSet(7,'Comic Sans MS');
            title('Interbeat Intervals (IBI) in samples')
            ylabel('IBI');
            set(h, 'MarkerSize', 15);
            xlim(xl);  
            grid on;
            
       %.............................................................
       % First Derivative
       %.............................................................
       xl = xlim;
            subplot(4,1,3)
            hold on;
            h = plot(fdt(MP), FD(MP), 'r.');
            AxisSet(7,'Comic Sans MS');
            title('First Derivative');
            ylabel('dx/dt');
            set(h, 'MarkerSize', 15);
            xlim(xl);
            grid on;
            
       
        %.............................................................
        % Second Derivative
        %.............................................................
        xl = xlim;
            subplot(4,1,4)
            hold on;
            h = plot(sdt(MP), SD(MP), 'r.');
            AxisSet(7,'Comic Sans MS');
            title('Second Derivative');
            xlabel('Samples');
            ylabel('d(dx/dt)/dt');
            set(h, 'MarkerSize', 15);
            xlim(xl);
            grid on;
            
    %-----------------------------------------------------------------
    % 2 on the keyboard: Rezoom
    %-----------------------------------------------------------------            
    elseif button == 50                              
        fprintf('\t\t=     Pausing for zoom. Zoom in then hit enter to continue   =\n');
        
        %.............................................................
        % Data and Peaks
        %.............................................................
        subplot(4,1,1)
            h = plot(k, x, 'b', k(MP), x(MP), 'r.');
            set(h, 'Markersize', 15);
            AxisSet(7,'Comic Sans MS');
            title('Manual Detection');
            ylabel('Magnitude');
            set(h, 'MarkerSize', 15);  
            axis(T) 
            grid on;
       
        zoom off;

       %.............................................................
       % Interbeat Intervals
       %.............................................................
       MP = sort(MP);
       P2P = diff(MP);
       dif = P2P;
       p2pi= (MP(1:length(MP)-1) + MP(2:length(MP)))/2;
       xl = xlim;
            subplot(4,1,2)
            h=plot(p2pi, dif);
            AxisSet(7,'Comic Sans MS');
            title('Interbeat Intervals (IBI) in samples')
            ylabel('IBI');
            set(h, 'MarkerSize', 15);
            xlim(xl); 
            grid on;
            
       %.............................................................
       % First Derivative
       %.............................................................
       xl = xlim;
            subplot(4,1,3)
            hold on;
            h = plot(fdt, FD, fdt(MP), FD(MP), 'r.');
            AxisSet(7,'Comic Sans MS');
            title('First Derivative');
            ylabel('dx/dt');
            set(h, 'MarkerSize', 15);
            xlim(xl);
            grid on;     
            
       %.............................................................
       % Second Derivative
       %.............................................................
        xl = xlim;
            subplot(4,1,4)
            hold on;
            h = plot(sdt, SD, sdt(MP), SD(MP), 'r.');
            title('Second Derivative');
            xlabel('Samples');
            ylabel('d(dx/dt)/dt');
            set(h, 'MarkerSize', 15);
            xlim(xl);
            grid on;
            pause; 
            
    %-----------------------------------------------------------------
    % 3 on the Keyboard: Remove Peak
    %-----------------------------------------------------------------              
    elseif button == 3
        
        %.............................................................
       % Data and Detection
       %.............................................................
        [mag, ind] = min(abs(mi-PI));
        pi = PI(ind);
        ri = find(MP~=pi);
        MP = MP(ri);
        MP = sort(MP);
        xl = xlim;
            subplot(4,1,1)
            hold off;
            h = plot(k, x, 'b', MP, x(MP), 'r.');
            title('Manual Detection');
            ylabel('Magnitude');
            set(h, 'MarkerSize', 15);  
            xlim(xl);
            grid on;
            
       %.............................................................
       % Interbeat Interval
       %.............................................................
       MP = sort(MP);
       P2P = [];
       P2P = diff(MP);
       dif = P2P;
       p2pi= (MP(1:length(MP)-1) + MP(2:length(MP)))/2;
       xl = xlim;
            subplot(4,1,2)
            hold off;
            h=plot(p2pi, dif);
            AxisSet(7,'Comic Sans MS');
            title('Interbeat Intervals (IBI) in samples')
            ylabel('IBI');
            set(h, 'MarkerSize', 15);
            xlim(xl);  
            hold on;
            grid on;
           
            
       %.............................................................
       % First Derivative
       %.............................................................
       xl = xlim;
            subplot(4,1,3)
            hold off;
            h = plot(fdt, FD, fdt(MP), FD(MP), 'r.');
            hold on;
            AxisSet(7,'Comic Sans MS');
            title('First Derivative');
            ylabel('dx/dt');
            set(h, 'MarkerSize', 15);
            xlim(xl);
            grid on;
                 
       
       %.............................................................
       % Second Derivative
       %.............................................................
        xl = xlim;
            subplot(4,1,4)
            hold off;
            h = plot(sdt, SD, sdt(MP), SD(MP), 'r.');
            hold on;
            title('Second Derivative');
            xlabel('Samples');
            ylabel('d(dx/dt)/dt');
            set(h, 'MarkerSize', 15);
            xlim(xl);
            grid on;
            hold on;
   
   %-----------------------------------------------------------------
   % 1 on the Keyboad: Replot
   %-----------------------------------------------------------------
   elseif button == 49                              
       MP = sort(MP);
       P2P = diff(MP);
       dif = P2P;
       p2pi= (MP(1:length(MP)-1) + MP(2:length(MP)))/2;
       xl = xlim;
       hold off
       subplot(4,1,2)
            h=plot(p2pi, dif);
            AxisSet(7,'Comic Sans MS');
            title('Interbeat Intervals (IBI) in samples')
            ylabel('IBI');
            set(h, 'MarkerSize', 15);
            xlim(xl);          
            
  %-----------------------------------------------------------------
  % 0 on the Keyboad: Close
  %-----------------------------------------------------------------                
  elseif button == 48                              
       
       mp  = MP(:);                                 
       ibi = diff(mp);
       close all;
       stop = 1;
       
   end;
    
end;

%=====================================================================
% Take care of outputs
%=====================================================================
if nargout==0,
    clear('mp', 'ibi');
end;  
