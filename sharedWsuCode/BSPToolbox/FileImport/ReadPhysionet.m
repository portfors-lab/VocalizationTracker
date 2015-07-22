function [Annotations] = ReadPhysionet(fileName)
%clear all;
%close all;
%fn = 'R:\PDAS\TBI\Patient007\Session001\Annotations.txt';

if ~exist(fileName,'file')
    error('Could not find file.');
end

dotIndex = find(fileName=='.');
record = fileName(1:dotIndex-1);
annotator = fileName(dotIndex+1:end);

executablePath = which('rdann.exe');
st = sprintf('!%s -r %s -a %s > TempPhysionetHeader.txt',executablePath,record,annotator);
eval(st);

fileIdentifier = fopen('TempPhysionetHeader.txt','r');
if fileIdentifier==-1
    error('Failed to open temporary physionet header text file.');
end

nAnnotationsIncrement = 1000;
Annotations = struct(...
    'Index'  ,zeros(nAnnotationsIncrement,1),...
    'Type'   ,zeros(nAnnotationsIncrement,1),...
    'SubType',zeros(nAnnotationsIncrement,1),...
    'Channel',zeros(nAnnotationsIncrement,1),...
    'Number' ,zeros(nAnnotationsIncrement,1),...
    'Comment',{cell(nAnnotationsIncrement,1)}...
    );
    
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
    
    Annotations.Index(cAnnotations) = str2num(line(13:21));
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

fileIdentifier = fclose(fileIdentifier);
if fileIdentifier==-1,
    error('Faled to close temporary physionet header text file.');
end

delete('TempPhysionetHeader.txt');

    