function [Data] = ReadPhysionetData(directory,Header)
%ReadPhysionetData: Reads data files from the Physionet
%
%   Data = ReadPhysionetData(directory,record,extension)
%   
%   directory     Directory where the annoation file is stored.
%   Header        Header structure returned by ReadPhysionetHeader.
%
%   Data          Structure containing vectors of the index, type, 
%                 subtype, channel, number, and comments for each 
%                 annotation.
%
%   This uses the Windows executable function rdann to create a 
%   temporary file that contains the annotations in a text file. This 
%   temporary file is then parsed and stored in the annotations 
%   structure and then deleted.
%   
%   Example: Read the annotations stored in the file named 
%   'mgh250.ari' of the MGH data base.
%
%      Annotations = ReadPhysionetAnnotations('R:\Physionet\MGH\Raw','mgh006.ari');
% 
%   A.L. Goldberger, L.A.N. Amaral, L. Glass, J.M. Hausdorff, 
%   P.C. Ivanov, R.G. Mark, J.E. Mietus, G.B. Moody, C. Peng, 
%   H.E. Stanley, "PhysioBank, PhysioToolkit, and PhysioNet: 
%   Components of a New Research Resource for Complex Physiologic 
%   Signals Circulation," 2000 , Vol. 101, pp. e215-e220.
%
%   Version 1.00 JM
%
%   See also ReadPhysionetData.

%====================================================================
% Error Checking
%====================================================================    
if nargin<2,
    help ReadPhysionetAnnotations;
    return;
    end;

if ~exist([directory '/' Header.RecordName '.dat'],'file')
    error('Could not find file.');
end

%====================================================================
% Author-Specified Parameters
%====================================================================    
tempFileName = 'TempPhysionetData.txt';

%====================================================================
% Preprocessing
%====================================================================    
record = Header.RecordName;
nChannels = Header.Channels;

%====================================================================
% Convert the Binary format to a Text File with rdann
%====================================================================
currentDirectory = pwd;
cd(directory);
executablePath = which('rdsamp.exe');
st = sprintf('!%s %s > %s',executablePath,record,tempFileName);
eval(st);
cd(currentDirectory);

%====================================================================
% Load the Temporary File and Delete It
%====================================================================
Data = load([directory '\' tempFileName]);
delete([directory '\' tempFileName]);

%====================================================================
% Post Processing
%====================================================================
Data = Data(:,2:end);

for c1=1:nChannels
    if Header.Gain(c1)~=0
        Data(:,c1) = (Data(:,c1)-Header.Baseline(c1))/Header.Gain(c1);
    else
        Data(:,c1) =  Data(:,c1)-Header.Baseline(c1);
    end
end


    