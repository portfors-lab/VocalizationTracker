function [X,ed,fs,F] = ReadAOFiles(fd);
%ReadAOFiles: Reads .mat files created by the Alpha Omega workstation
%
%   [X,fs,ed,fn] = ReadAOFiles(fd)
%   
%   fd   Name of directory containing MER recordings.
%
%   X    Cell array of vectors containing the MER recordings.
%   fs   Vector of sample rates for each recording.
%   ed   Electrode depth for each recording.
%   F    Cell array of file names for each recording.
%
%   This file reads .mat files created by the Alpha Omega 
%   surgical workstation. Once the files are read in, I recommend
%   that the three structures be stored in a single .mat file. 
%   This file will be much smaller and quicker to read than the 
%   individual .mat files, which contain a lot of superfluous 
%   information.
%   
%   Example: Read all of the data for STN129R.
%
%      pc = 'STN129R';
%      t0 = clock;
%      if exist([pc '.mat'],'file'),
%         load([pc '.mat']);
%         fprintf('Loaded from local .mat file in %d seconds.\n',round(etime(clock,t0)));
%      else
%         [X,ed,fs,F] = ReadAOFiles(['R:\Microelectrode\Burchiel\STN\' pc]);
%         save(pc,'X','ed','fs','F');
%         fprintf('Loaded from AO .mat files in %d seconds.\n',round(etime(clock,t0)));
%         end;
% 
%   J. H. Falkenberg, J. McNames, K. J. Burchiel, "Statistical 
%   methods of analysis and visualization of extracellular 
%   microelectrode recordings," Annual International Conference of 
%   the IEEE Engineering in Medicine and Biology - Proceedings, 
%   Cancun, Mexico, pp. 2515-2518, 17-21 September 2003.
%
%   Version 1.00 JM
%
%   See also NonparametricSpectrogram.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help ReadAOFiles;
    return;
    end;
    
fl = dir([fd '/*.mat']);                                   % File list

if length(fl)==0,
    error('No data files found in current directory.');
    end;   
    
%====================================================================
% Preprocessing
%====================================================================
nf = length(fl);                                           % Number of files
cf = 0;                                                    % Number of valid files read

%====================================================================
% Memory Allocation
%====================================================================
X  = cell(nf,1);
fs = zeros(nf,1);
ed = zeros(nf,1);
F  = cell(nf,1);

%====================================================================
% Load the Data
%====================================================================
for c1=1:nf,
    fn = fl(c1).name;                                      % File name
    ao = load([fd '/' fn]);                                % Load the AO file structure
    
    if ~isfield(ao,'CElectrode1_KHz'),
        warning(sprintf('File %s appears to be empty. Skipping.',fn));
        continue;
        end;
        
    if length(ao.CElectrode1)==0,
        warning(sprintf('Apparently no MER data recorded. Skipping.',fn));
        continue;
        end;
        
    cf     = cf + 1;    
    X{cf}  = ao.CElectrode1;                               % MER data
    ed(cf) = str2double(fn(1:5))/1e3-11;                   % Electrode depth
    fs(cf) = round(ao.CElectrode1_KHz*1e3);                % Sample rate
    F{cf}  = fn;                                           % File name
    drawnow;
    end;  
    
%====================================================================
% Postprocessing
%====================================================================
X  = X(1:cf);
ed = ed(1:cf);
fs = fs(1:cf);
F  = F(1:cf);

%====================================================================
% Output Arguments
%====================================================================
if nargout==0,
    clear('X','ed','fs','F');
    end;