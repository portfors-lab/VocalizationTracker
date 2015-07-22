function [pinkNoise transient] = Pinkify(y)
% Generate pink noise from white, gaussian noise.

nSamples = length(y);
B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
A = [1 -2.494956002   2.017265875  -0.522189400];
pinkNoise = filter(B,A,y);    % Apply 1/F roll-off to PSD
transient = round(log(1000)/(1-max(abs(roots(A)))));
% pinkNoise(1) = 0;
% pinkNoise(end) = 0;


% nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est. = 1430

