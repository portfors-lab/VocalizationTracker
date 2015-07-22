function [x_p] = Percentile(x,P,arg1);
% Percentile(x,P,[N])
% Estimate the P'th percentile value for x
% P should be a value between 0 and 1.  
% If P is a vector the function returns a
% vector of the same size with the corresponding
% percentiles.
N = length(x);
if exist('arg1'),
   N = arg1;
	end;
   
x_p = zeros(size(P));
sx = sort(x(1:N));
for cnt = 1:length(P),
	ti       = ceil(P(cnt)*N);
	if ti>1,
		x_p(cnt) = mean([sx(ti-1) sx(ti)]);
	else
		x_p(cnt) = sx(ti);
		end;
	end;
