function VocalizationTracker

welcomepicture = imread('batmouse.jpg');
fig = figure('units', 'normalized', 'position', [0.25 0.15 0.6 0.7], 'menubar', 'none');


 image(welcomepicture)
 axis off
axis image

pos = get(gca, 'position');
set(gca, 'position', [pos(1)+0.2 pos(2) pos(3) pos(4)])

uicontrol(fig, 'style', 'text', 'units', 'normalized', 'position', [0.05 0.8 0.4 0.1], 'string', 'Vocalization Tracker', 'backgroundcolor', [0.8 0.8 0.8], 'fontsize', 18)

uicontrol(fig, 'style', 'pushbutton', 'units', 'normalized', 'position', [0.1 0.6 0.25 0.05], 'string', 'New Spectrogram', 'callback', @newSpec);
uicontrol(fig, 'style', 'pushbutton', 'units', 'normalized', 'position', [0.1 0.4 0.25 0.05], 'string', 'Load Previous Session', 'callback', @LoadPreviousSession);


function LoadPreviousSession(ho, ed)

[f p] = uigetfile('*.mat');
h=[];
userInputs=[];
filePath=[];
plotdata = [];
load([p f]);

MainSpecFigure(filePath, userInputs, h, plotdata);
close(fig)
end

function newSpec(ho, ed)

MainSpecFigure
close(fig)
end 

end