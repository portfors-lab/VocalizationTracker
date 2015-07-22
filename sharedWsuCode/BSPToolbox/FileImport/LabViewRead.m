function [X,t] = LVRead(fn,nxa,pf);
%LVRead: Reads LabView data.
%
%   [X,t] = LVRead(fn,nx,pf);
%
%   fn   File name (string).
%   nx   Number of samples to read. Default = all.
%   pf   Plot flag. 
%
%   X    All numeric signal data in the file.
%   t    Sample times of each signal.
%
%   Reads LabView data files provided by Prasad's lab. This could be
%   significantly improved to automatically read more information 
%   about the files such as the units, time collected, and sample
%   rate. 
% 
%   This was developed based on data acquired Dr. Prasad's lab. 
% 
%   If PlotFlag is included, will generate plot of each lead.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help LVRead;
    return;
    end;

nx = inf; 
if exist('nxa','var') & ~isempty(nxa),
    nx = nxa;
    end;        
    
%====================================================================
% Author-Specified Parameters
%====================================================================  
as = 100000;                                               % Allocation size for new data    
    
%====================================================================
% Process function arguments
%====================================================================  
pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Open File
%====================================================================  
fid = fopen(fn,'r');
if fid==-1,
    error('Could not open file %s.\n',fn);
end

%====================================================================
% Skip to Data Fields
%====================================================================  
fl = fgetl(fid);  
ln = 1;                                                    % Line number
while isempty(fl) | fl~=-1,
    d = str2num(fl);
    fl = fgetl(fid);          
    ln = ln + 1;                                           
    if size(d,1)==1 | size(d,2)>=2,                        % Look for a line with at least two numbers
        break;
        end;
    end;

if fl==-1,
    fclose(fid);
    error('Could not find any valid data in the file.');
    end;

if rem(d,2)~=0,
    fclose(fid);
    error('File contains odd number of columns.');
    end;

%====================================================================
% Allocate Memory
%====================================================================  
nc = length(d);                                            % Number of columns
D  = zeros(as,nc);

nd = 1;
D(nd,:) = d(1:end);

%====================================================================
% Main Loop
%====================================================================  
while fl~=-1, 
    d = str2num(fl);
    fl = fgetl(fid);
    ln = ln + 1;
    if size(d,1)~=1 | size(d,2)~=nc,
        fclose(fid);
        error('Encountered invalid data after start of data.');
        end;
    if nd+1>size(D,1),
        Dt = D;
        D  = zeros(size(D,1)+as,nc);
        D(1:nd,:) = Dt;
        fprintf('Expanding array size...\n');
        drawnow;
        end;
    nd = nd + 1;
    D(nd,:) = d;   
    if rem(nd,10000)==0,
        fprintf('Working on %d\n',nd);
        drawnow;
        end;
    if nd==nx,
        break;
        end;
    end;

fclose(fid);

if size(D,2)==0,
    fprintf('What happened?\n');
    pause;
    end;

%====================================================================
% Postprocessing
%====================================================================  
D = D(1:nd,:);                                             % Eliminate all zero padding for memory allocation
t = D(:,1);                                                % Get the time index
xi = [];
for c1=2:nc,                                               % Determine which columns are data, not time indices
    if ~all(D(:,c1)==t),
        xi = [xi;c1];
        end;
    end;
X = D(:,xi);
    
%====================================================================
% Plot Signal That Was Read
%====================================================================
if pf==1,
    figure;
    FigureSet;
    plot(t,X);
    xlabel('Time (sec)');
    ylabel('Signals');
    title('LVM Signals');
    box off;
    AxisSet;
    end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('X','t');
    end;
