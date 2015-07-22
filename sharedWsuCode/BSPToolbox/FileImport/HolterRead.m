function [y,fs] = HolterRead(fn,lda,npa,pf);
%HolterRead: Reads Holter ECG data.
%
%   [y,fs] = HolterRead(fn,ld,np,pf);
%
%   fn   File name (string).
%   ld   Lead (integer, 1-3). Default = 1.
%   np   Number of points to read. Default = inf.
%        If a vector, specifies first and last point to read.
%   pf   Plot flag. 
%
%   y    Requested signal segment of the electrocardiogram (scaled).
%   fs   Sample rate (Hz).
%
%   Reads the Holter waveform data for the filename specified. 
%   Returns a matrix with three columns where each column contains 
%   the QRS waveform from each of the three leads. The sample rate is 
%   fixed at 128 Hz.
% 
%   This was developed based on data acquired from a Holter monitor
%   and some careful reverse-engineering of the file format. It may 
%   not work with all data acquired from Holter monitors.
%
%   If the size is specified, this only reads in the first Size 
%   samples. Otherwise, it reads the entire file.
% 
%   If PlotFlag is included, will generate plot of each lead.

%====================================================================
% Error Checking
%====================================================================    
if nargin<2,
    help HolterRead;
    return;
    end;

%====================================================================
% Process function arguments
%====================================================================
ld = 1;
if exist('lda') & ~isempty(lda),
    ld = lda;
    end;

np = [1 inf];
if exist('npa') & ~isempty(npa),
    if length(npa)==1,
        np = [1 max(1,npa)];
    else
        np = npa;
        end;
    end;
   
pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Author-Specified Parameters
%====================================================================  
ns = np(1)-1;       % No. of points to skip forward
nr = np(2)-np(1)+1; % No. of points requested
fs = 128;           % Sampling rate in HZ 
ng = 115200;        % Number of test points in the file that are garbage

%====================================================================
% Main Routine
%====================================================================  
fid = fopen(fn,'r');
if fid==-1,
    fprintf('Error: Could not open file %s.\n',fn);
    return;
    end
        
sf = fseek(fid,ng,'bof');    % Skip over beginnig configuration segment   
if sf==-1,
    fprintf('Error: Could reach start of data in file %s.\n',fn);
    return;
    end    
    
sf = fseek(fid,ld-1,'cof');  % Skip forward to specified lead     
if sf==-1,
    fprintf('Error: Could reach start of lead data in file %s.\n',fn);
    return;
    end

sf = fseek(fid,3*ns,'cof');  % Skip forward to user-specified start of segment    
if sf==-1,
    fprintf('Error: Could reach start of lead data in file %s.\n',fn);
    return;
    end    
    
    
[y,cr] = fread(fid,nr,'char',2); % Read specified segment of data

fclose(fid);

if nr~=cr & nr<inf,
    fprintf('Warning: Only %d points read of %d requested.\n',cr,nr);
    end;
    
%====================================================================
% Postprocessing
%====================================================================  
y = y(:); % Convert to column vector    
    
%====================================================================
% Plot Signal That Was Read
%====================================================================
if pf==1,
    ny = length(y);
    k  = 1:ny;
    t  = (k-1)/fs;
    figure;
    FigureSet;
    plot(t,y(:,1));
    xlabel('Time (sec)');
    ylabel('ECG Signal (Scaled)');
    title('Holter ECG Waveform');
    AxisSet;
    zoom on;
    end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('y','fs');
    end;
