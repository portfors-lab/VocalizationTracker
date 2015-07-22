function result = ismac
%ISMAC True for the Mac OS X version of MATLAB.
%   ISMAC returns 1 for MAC (Macintosh) versions of MATLAB and 0 otherwise.
%
%   See also COMPUTER, ISUNIX.

%   Copyright 1984-2006 The MathWorks, Inc. 
%   $Revision: 1.1.1.1 $  $Date: 2011-04-19 17:36:44 $

result = strncmp(computer,'MAC',3);