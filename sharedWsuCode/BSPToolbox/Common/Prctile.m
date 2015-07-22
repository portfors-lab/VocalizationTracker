function [y] = Prctile(x,p,arg1);
%Prctile: Calculate percentile(s) of a signal x
%
%   [y] = Prctile(x,P,arg1);
%
%   x    Input vector.
%   p    Vector of percentiles (%). 
%
%   y    Percentiles of input vector x. 
%
%   Serves as a replacement of prctile for people who don't have
%   the statistics toolbox.
%
%   Example: Calculate the 5th and 95th percentiles of a random
%   vector.
%
%      x = rand(500,1);
%      p = Prctile(x,[5 95]);
%
%   M. Hayes, Statistical Digital Signal Processing and Modeling. 
%   New York: John Wiley & Sons, 1996, pp. 415-420.
%
%   Version 0.00.01 JM
%
%   See also Common.

n  = length(x); 
y  = zeros(size(p));
xs = sort(x);
for c1 = 1:length(p),
    pt    = max(p(c1)/100,0);
    pt    = min(pt,100);
	id    = pt*(n-1) + 1;
    lid   = floor(id);
    uid   = ceil(id);
    if uid==lid,
        y(c1) = xs(id);
    else
        y(c1) = (id-lid)*(xs(uid)-xs(lid))/(uid-lid) + xs(lid);
        end;
	end;
