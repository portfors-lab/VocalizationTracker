function [pis,fom,tmp] = TemplateDetect(x,fs,ipi,tdr,ftya,pfa);
%TemplateDetect: Uses the template matching method to locate peaks
%
%   [pis,fom,tmp] = TemplateDetect(x,fs,ipi,tdr,fty,pf)
%
%   x     Signal in which events are to be detected.
%   fs    Sample rate (Hz).
%   ipi   Initial event indices to estimate the template.
%   tdr   2-dimensional vector of the number of samples before
%         and after event indices that represent the template 
%         duration
%   ftya  Figure of merit type argument: 1=MSE (default), 
%         2=Cross-correlation.
%   pfa   Plot flag argument: 0=none (default), 1=screen
%
%   pis   Sorted indices of locations where events were detected.
%   fom   Figure of merit
%   tmp   Template
%
%   To be written.
%
%   Version 1.00 SK
%
%   See also DetectMaxima, DetectMinima, and PowerPeaks.

%====================================================================
% Error Checking
%====================================================================    
if nargin<4,
    help TemplateDetect;
    return;
    end;
   
%====================================================================
% Process Function Arguments
%====================================================================
fty = 2;
if exist('ftya') & ~isempty(ftya),
    fty = ftya;
    end;   

pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
frmin=10;
frmax=250;
   
%====================================================================
% Preprocessing
%====================================================================
nmina = ceil(frmin*(length(x)-1)/fs);                      % Minimum number of spikes
nmaxa = ceil (frmax*(length(x)-1)/fs);                     % Maximum number of spikes

x   = x(:);                                                % Convert signal to column vector, if not already
k   = -tdr(1):tdr(2);                                      % Row vector of offset indices
%dk  = round(-2*tdr(1)):round(2*tdr(2));                    % Row vector of offset indices just for displaying
ntp = length(k);                                           % Number of template points
ipi = ipi(:);                                              % Convert peak indices to column vector, if not already
nip = length(ipi);                                         % Number of initial points
nx  = length(x);                                           % Signal length (samples)
%ndtp= length(dk);

%====================================================================
% Estimate the Template
%====================================================================
SI  = ipi*ones(1,ntp) + ones(nip,1)*k;                     % Matrix of signal indices : one row for one spike
vri = find(all(SI>0 & SI<=nx,2));                          % Identify the valid row indices (those than don't exceed the signal limits)
SI  = SI(vri,:);                                           % Strip out invalid rows
%DSI = ipi(1:end-1,1)*ones(1,ndtp) + ones(nip-1,1)*dk;      % Matrix of signal indices just for displaying
tmp = median(x(SI));                                       % Estimated template vector(row)

%====================================================================
% Calculate Figure of Merit
%====================================================================
switch fty,
    case 1,                                                % Mean Squared Error
        fom = TemplateError(x,fs,tmp);                     % Calculated template error, i.e., the Euclidean distance
        emx = max(fom);                                    % Calculate the maximum template error
        fom = [emx*ones(tdr(1),1);fom;emx*ones(tdr(2),1)]; % Pad with maximum error to deal with edge effects
	    pi  = DetectMinima(fom,[],0);                      % Only count minima as candidate peaks
        [jnk,id] = sort(fom(pi),1,'ascend');               % Sort the best FOM of points in ascending order        
    case 2,                                                % Matched filter
		b   = tmp(end:-1:1);                               % Numerator filter coefficients(flip the order of the template)
        a   = 1;                                           % Denominator filter coefficients
		fom = filter(b,a,x);                               % Use a filter function to calculate the crosscorrelation
        ns  = tdr(2);                                      % Number of points to shift         
        fom = [fom(ns+1:end,:);zeros(ns,1)];               % Shift FOM so that it is alligned with the signal      
        pi  = DetectMaxima(fom,[],0);                      % Only count maxima as candidate peaks
        [jnk,id] = sort(fom(pi),1,'descend');              % Sort the best FOM of points in descending order
    end;
pis = pi(id);

%====================================================================
% Generate Plot to Show Analysis
%====================================================================    
% Do not generate figures for now
if 0,   
	switch fty,
        case 1,
            pis = pis(pis>0 & pis<=nx);                    % Strip out invalid peak indices that exceed the signal range
            k = 1:length(x);
            t = (k-1)/fs;
            thmin= fom(pis(nmina));
            th   = thmin;
            ipi  = pis(find(fom(pis)<thmin));              % The largest peak indices
            
            figure;
            FigureSet;
                h = plot(t(1:length(fom)),fom,'k',t(ipi),fom(ipi),'ko');
                set(h(1),'Color',0.7*[1 1 1]);
                set(h(2),'MarkerFaceColor',0.2*[1 1 1]);
                set(h(2),'MarkerEdgeColor',1*[1 1 1]);
                set(h(2),'MarkerSize',7);
                hold on;
                    g=plot(t,th*ones(length(t),1),'k-.');
                    set(g,'LineWidth',2);
                hold off;
                box off;
                xlabel('Time (s)');
                ylabel('FOM (scaled)');
                xlim([0 0.5]);
                ylim([min(1.025*fom(1:ceil(fs/2))) max(1.025*fom(1:ceil(fs/2)))]);
                FormatTicks('%1.2f','%1.2f');
                AxisSet(8,'Times New Roman');
                set(gca,'Layer','top');                
        case 2,
            pis = pis(pis>0 & pis<=nx);                    % Strip out invalid peak indices that exceed the signal range
            k = 1:length(x);
            t = (k-1)/fs;
            thmax= fom(pis(nmina));
            th   = thmax;
            ipi  = pis(find(fom(pis)>thmax));               % The largest peak indices
      
            figure;
            FigureSet;
                h = plot(t,fom,'k',t(ipi),fom(ipi),'ko');
                set(h(1),'Color',0.7*[1 1 1]);
                set(h(2),'MarkerFaceColor',0.2*[1 1 1]);
                set(h(2),'MarkerEdgeColor',1*[1 1 1]);
                set(h(2),'MarkerSize',7);
                hold on;
                    g=plot(t,th*ones(length(t),1),'k-.');
                    set(g,'LineWidth',2);
                hold off;
                box off;
                xlabel('Time (s)');
                ylabel('FOM (scaled)');
                xlim([0 0.5]);
                legend('Figure of merit','Peaks','Threshold');
                FormatTicks('%1.2f','%1.2f');
                AxisSet(8,'Times New Roman');
                set(gca,'Layer','top');
    end;
end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('pi');
    end;

