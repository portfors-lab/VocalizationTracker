function [c] = DistinctColors(nColors);
%DistinctColors: Generates a set of distinct colors
%
%   [c] = DistinctColors();
%
%   c    Matrix containing distinct colors (rgb pairs in columns).
%
%   Generates colors that are visually distinctive. Is useful for
%   representing categorical data in two-dimensional plots by
%   color. 
%   
%   Example: Generate the spectrogram of an intracranial pressure
%   signal using a Blackman-Harris window that is 45 s in duration.
%
%      x = randn(50,5);
%      y = randn(size(x));
%      c = DistinctColors;
%      h = plot(x,y,'.');
%      for c1=1:size(x,2),
%         set(h(c1),'Color',c(c1,:));
%         end;
%
%   C. Ware, "Information Visualization: Percetpion for Design," 2nd edition,
%   Morgan Kaufmann, 2004.
%
%   Version 1.00 JM
%
%   See also ColorSpiral. 

c = ...
[1.00 0.00 0.00;... % Red
 0.00 0.00 1.00;... % Blue
 0.00 1.00 0.00;... % Green
 1.00 1.00 0.00;... % Yellow
 0.00 0.00 0.00;... % Black
 1.00 0.00 1.00;... % Violet
 0.00 1.00 1.00;... % Cyan
 1.00 0.50 0.25;... % Orange 
 0.50 0.50 0.25;... % Olive 
 0.50 0.10 0.50;... % Dark Violet
 0.10 0.50 0.50;... % Dark cyan
 0.70 0.20 0.75;... % Magenta
 0.75 0.33 0.10;... % Brown
 0.11 0.83 0.83;... % Sky blue
 0.50 0.50 0.50;... % Gray
 1.00 0.50 0.50;... % Salmon
 0.50 0.10 1.00;... % Purple
 0.25 0.50 0.50;... % Blue-green
 0.50 0.10 0.00;... % Dark Red
 0.75 0.91 0.06;... % Lime
 0.87 0.61 0.06;... % Pumpkin
 0.55 0.22 0.71;... % Grape
 0.84 0.82 0.71;... % Tan
];

if exist('nColors','var')
    c = c(1:nColors,:);
end
