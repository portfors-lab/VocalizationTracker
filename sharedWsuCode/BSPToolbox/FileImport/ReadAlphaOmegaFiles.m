function [X] = ReadAlphaOmegaFiles(fd,channelNameArgument,dfa)
%ReadAOFiles: Reads .mat files created by the Alpha Omega workstation
%
%   [X,fs,ed,fn] = ReadAlphaOmegaFiles(fd,sn)
%   
%   fd     Name of directory containing MER recordings.
%   sn     Cell array of signal names. Default = {'MER1'};
%   tmin   Minimum duration of the recording. Default = 0.
%   df     Display flag. Default=1.
%
%   X      Structure array containing the data.
%   nf     Number of files found in the specified directory.
%
%   This file reads .mat files created by the Alpha Omega 
%   surgical workstation. Each structure includes a field name
%   for the requested signals (e.g., MER1, MER2, EMG1, EMG2,
%   LFP1, LFP2), the filename, and the electrode depth. Each
%   signal is also contained in a structe that includes
%
%      x          A vector containing the signal.
%      fs         The sample rate (Hz).
%      Signal     Signal name.
%      TimeBegin  The time at which the recording began.
%      TimeEnd    The time at which the recording ended. 
%
%   Note that this function uses a different naming scheme for the
%   fields than the internal format of the Alpha Omega .mat files.
%   Here is the name conversion mapping.
% 
%      MER1--CElectrode1   LFP1--LFP1   EMG1--EMG1
%
%   Once the files are read in, I recommend
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
    help ReadAlphaOmegaFiles;
    return;
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================    
channelName = {'MER1'};
if exist('channelNameArgument','var') && ~isempty(channelNameArgument),
    channelName = channelNameArgument;
    end;         

displayFlag = 1;
if exist('dfa','var') && ~isempty(dfa), %#ok 
    displayFlag = dfa;    
    end;         
    
%====================================================================
% Preprocessing
%====================================================================    
fl = dir([fd '/*.mat']);                                   % File list

if length(fl)==0,
    error('No data files found in current directory.');
    end;   
    
%====================================================================
% Preprocessing
%====================================================================
nFiles      = length(fl);                                  % Number of files
cFilesValid = 0;                                           % Number of valid files read

%====================================================================
% Memory Allocation
%====================================================================
nChannels       = length(channelName);                     % Number of channels
%electrodeDepths = zeros(n,1);
%validFiles      = zeros(nf,1);

%cs = struct('x',[],'fs',nan,'Signal','','TimeBegin',nan,'TimeEnd',nan); % Channel structure
%X  = repmat(struct('FileName','','ElectrodeDepth',[],'Channel',repmat(cs,nc,1)),nf,1);  % Recording structure

X = struct(...
    'PatientCode','',...                                   % Patient Code
    'TrajectoryName','',...                                % Trajectory Name
    'NoFiles',nFiles,...                                   % No. files
    'NoChannels',nChannels,...                             % No. channels
    'ElectrodeDepth',zeros(nFiles,1),...                   % Electrode depths
    'FileName',[],...                                      % File Names
    'Channel',repmat(struct(...
        'Name','',...                                      % Name of the signal recorded by this channel
        'Description','',...                               % Signal Description
        'SampleRate',nan,...                               % Sample rate (Hz)
        'Recording',repmat(struct(...
            'Signal',[],...                                % Signal vector
            'TimeBegin',nan,...                            % Time of start of recording (s)
            'TimeEnd',nan...                               % Time of end of recording (s)
            ),nFiles,1)...                                 % Cell array of signal vectors
        ),nChannels,1)...
    );
X.FileName{nFiles} = '';                                   % Allocate memory for file names

%====================================================================
% Load the Data
%====================================================================
t0 = clock;
for c1=1:nFiles,
    fileName = fl(c1).name;                                % File name
    fileStructure = load([fd '/' fileName]);               % Load the AO file structure
    
    X.FileName{c1}       = fileName;
    X.ElectrodeDepth(c1) = str2num(fileName(1:find(fileName=='.')-1));
    
    for c2=1:length(channelName),                           % Loop over the requested signal names
        if strcmpi(channelName{c2}(1:3),'MER'),
            alphaOmegaField = sprintf('CElectrode%c',channelName{c2}(end));
        else
            alphaOmegaField = sprintf('C%s',channelName{c2});
            end;
            
        if ~isfield(fileStructure,alphaOmegaField),
            continue;
            end;
            
        x   = getfield(fileStructure,sprintf('%s'          ,alphaOmegaField));
        fs  = 1e3*getfield(fileStructure,sprintf('%s_KHz',alphaOmegaField));
        fso = 1e3*getfield(fileStructure,sprintf('%s_KHz_Orig',alphaOmegaField));      
        
        if c1>1 & fs~=X.Channel(c2).SampleRate,
            warning(sprintf('Sample rate in file %s differs from previous files.',fileName));
            end;
        if fs~=fso,
            warning(sprintf('Original and actual sample rates were different for %s in file %s.',alphaOmegaField,fileName));
            end;
            
        X.Channel(c2).Name                    = channelName{c2};
        X.Channel(c2).SampleRate              = fs;
        X.Channel(c2).Recording(c1).Signal    = x;
        X.Channel(c2).Recording(c1).TimeBegin = getfield(fileStructure,sprintf('%s_TimeBegin',alphaOmegaField));
        X.Channel(c2).Recording(c1).TimeEnd   = getfield(fileStructure,sprintf('%s_TimeEnd'  ,alphaOmegaField));        
        end;
    if displayFlag,
        t1 = clock;                                        % Current time
        te = etime(t1,t0);                                 % Time elapsed
        tr = (te/c1)*(nFiles-c1);                          % Time remaining
        fprintf('(%3d of %3d) Elapsed:%4.1f s   Remaining:%4.1f s   %s\n',c1,nFiles,te,tr,fileName);
        end;
    drawnow;
    end;  
    
%====================================================================
% Postprocessing
%====================================================================


%====================================================================
% Output Arguments
%====================================================================
if nargout==0,
    clear('X');
    end;