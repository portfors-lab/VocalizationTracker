function [s] = IsScript(fn)

warning('IsScript is obsolete. Please remove from your code.');
s=0;

return;

%IsScript: Returns 1 if a .m file is a script, returns 0 otherwise
%
%   s = IsScript(fn)
%
%   fn   Name of function (text).       
%
%   s    Flag: 1 if is a script, 0 otherwise   
%
%   This is primarily designed to allow .m files to determine 
%   automatically whether they are being called as a script or not. 
%   This operates by scanning the first line in the file to determine
%   whether it begins with the word 'function' or not.
%
%   Example: Determine whether the file named 'Spectrogram.m' is 
%   a script or not.
%
%      s = IsScript('Spectrogram');
%
%   Version  1.00 JM
%
%   Desmond J. Higham and Nicholas J. Higham, "MATLAB Guide," 2000.
%
%   See also mfilename, which, and strncmpi.

%====================================================================
% Process function arguments
%====================================================================
if nargin~=1,
    help IsScript;
    return;
    end;

%====================================================================
% Main Code
%====================================================================  
fp = which(fn); % Obtain complete path to file name
fid = fopen(fp,'r');
if fid==-1,
    fprintf('ERROR: IsScript could not open specified file.');
    return;
    end;

fl = fgetl(fid); % Get the first line of the file
fclose(fid);

if strncmp('P-file',fl(2:7),6), % Assume .p files are not scripts
    s = 0;
    return;
    end;

s = ~strncmpi('function',fl,8);