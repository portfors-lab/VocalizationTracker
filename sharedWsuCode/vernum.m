function vernum = vernum();
% vernum  - returns the current version of Matlab as a number
%
% 30 June 2006 - Created
% by David S. Ciochetto and Audrey Barnett, Dalhousie Oceanography
% 30 June 2006 (Last modified) by David S. Ciochetto
%
% This function gets your current version of matlab and converts up through
% the first decimal point to a numeric value. It is intended to be used to
% set off cases where a fix is needed for old code for a script that you
% develop. This way the user does not have to edit the code if they are
% using an older version.
%
% Example
% if vernum < 7.1 % The below works for Matlab versions earlier than 7.1
%   % Add code
% else % works for newer code
%   % Add fancy code here
% end
%
% NOTES:
%
% USES:
%
% SEE ALSO:
%

% Select a method, other methods retained as notes
do = 'a';

% Get the version
tmp = version;

% Convert it to a number
if do == 'd'  % Dave's Version (Pascal minded)
  ii = 1;
  while ~isempty(str2num(tmp(1:ii)))
     vernum = str2num(tmp(1:ii));
     ii = ii + 1;
  end
elseif do == 'a'  % Audrey's version (sans loop)
  x=strfind(tmp,'.');
  vernum=str2num(tmp(1:x(2)-1));
end
