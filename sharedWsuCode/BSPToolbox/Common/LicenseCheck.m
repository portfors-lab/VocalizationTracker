function [] = LicenseCheck
%LicenseCheck: Check which users are using each of the toolboxes
%
%   [] = LicenseCheck
%
%   A simple utility to check who is currently using MATLAB. 
%
%   Example: check who is using MATLAB.
%
%      LicenseCheck
%
%   Desmond J. Higham and Nicholas J. Higham, "MATLAB Guide," 2000.
%
%   Version 1.00 JM
%
%   See also AxisSet and FormatTicks.

st = sprintf('!lmutil lmstat -a -c "%s\\bin\\win32\\license.dat"',matlabroot);
eval(st);