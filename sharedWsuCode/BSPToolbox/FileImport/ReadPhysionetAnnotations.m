function [Annotations] = ReadPhysionetAnnotations(directory,record,extension)
%ReadPhysionetAnnotations: Reads annotations files from the Physionet
%
%   Annotations = ReadPhysionetAnnotations(directory,record,extension)
%   
%   directory     Directory where the annoation file is stored.
%   record        Name of the record.
%   extension     3 character file extension
%
%   Annotations   Structure containing vectors of the index, type, 
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
%      Annotations = ReadPhysionetAnnotations('R:\Physionet\MGH\Raw','mgh006','ari');
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
if nargin<2
    help ReadPhysionetAnnotations;
    return;
end

if ~exist([directory '/' record '.' extension],'file')
    error('Could not find file.');
end

%====================================================================
% Author-Specified Parameters
%====================================================================    
tempFileName = 'TempPhysionetAnnotations.txt';

%====================================================================
% Convert the Binary format to a Text File with rdann
%====================================================================
currentDirectory = pwd;
cd(directory);
executablePath = which('rdann.exe');
st = sprintf('!%s -r %s -a %s > %s',executablePath,record,extension,tempFileName);
eval(st);
cd(currentDirectory);

%====================================================================
% Open the text file for reading
%====================================================================
fileIdentifier = fopen([directory '\' tempFileName],'r');
if fileIdentifier==-1
    error('Failed to open temporary physionet annoatations text file.');
end

%====================================================================
% Allocate memory for the annotations
%====================================================================
nAnnotationsIncrement = 1000;
Annotations = struct(...
    'Index'  ,zeros(nAnnotationsIncrement,1),...
    'Type'   ,zeros(nAnnotationsIncrement,1),...
    'SubType',zeros(nAnnotationsIncrement,1),...
    'Channel',zeros(nAnnotationsIncrement,1),...
    'Number' ,zeros(nAnnotationsIncrement,1),...
    'Comment',{cell(nAnnotationsIncrement,1)}...
    );
    
%====================================================================
% Main loop
%====================================================================    
line = fgetl(fileIdentifier);
cAnnotations = 0;
nBufferAnnotations = nAnnotationsIncrement;
while line~=-1,
    cAnnotations = cAnnotations + 1;
    if cAnnotations > nBufferAnnotations,
       nBufferAnnotations = nBufferAnnotations + nAnnotationsIncrement;
       AnnotationsTemp = struct(...
            'Index'  ,zeros(nBufferAnnotations,1),...
            'Type'   ,zeros(nBufferAnnotations,1),...
            'SubType',zeros(nBufferAnnotations,1),...
            'Channel',zeros(nBufferAnnotations,1),...
            'Number' ,zeros(nBufferAnnotations,1),...
            'Comment',{cell(nBufferAnnotations,1)}...
            );
        AnnotationsTemp.Index  (1:cAnnotations-1) = Annotations.Index;
        AnnotationsTemp.Type   (1:cAnnotations-1) = Annotations.Type;
        AnnotationsTemp.SubType(1:cAnnotations-1) = Annotations.SubType;
        AnnotationsTemp.Channel(1:cAnnotations-1) = Annotations.Channel;
        AnnotationsTemp.Number (1:cAnnotations-1) = Annotations.Number;
        AnnotationsTemp.Comment(1:cAnnotations-1) = Annotations.Comment;
        Annotations = AnnotationsTemp;
    end
    
    Annotations.Index(cAnnotations) = str2double(line(13:21));
    Annotations.Type (cAnnotations) = line(27);
    Annotations.SubType(cAnnotations) = line(32);
    Annotations.Channel(cAnnotations) = line(37);
    Annotations.Number (cAnnotations) = line(42);
    if length(line)>43
        Annotations.Comment{cAnnotations} = line(44:end);
    end
    line = fgetl(fileIdentifier);
end

Annotations.Index   = Annotations.Index  (1:cAnnotations);
Annotations.Type    = Annotations.Type   (1:cAnnotations);
Annotations.SubType = Annotations.SubType(1:cAnnotations);
Annotations.Channel = Annotations.Channel(1:cAnnotations);
Annotations.Number  = Annotations.Number (1:cAnnotations);
Annotations.Comment = Annotations.Comment(1:cAnnotations);

%====================================================================
% Close and delete the temporary text file
%====================================================================
fileIdentifier = fclose(fileIdentifier);
if fileIdentifier==-1,
    error('Faled to close temporary physionet annotations text file.');
end

delete([directory '\' tempFileName]);

    