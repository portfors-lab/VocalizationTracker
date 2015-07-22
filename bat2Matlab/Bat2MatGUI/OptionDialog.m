function [stimPath, colormap, invertColor, colorRange, txtPath] = OptionDialog(stimPath, colormap, invertColor, colorRange, txtPath)
% function [stimPath, colormap, colorRange] = optionDialog(stimPath, colormap, colorRange, colorChoice)
%
%INPUTS
%
%stimPath          folder where stimulus files for the data reside
%colormap          string of the desired colorMap to use
%colorRange        manual selection of range of colorbar. default = []
%colorChoice       indicates which radio button to have selected for
%                  colorRange

% Author Amy Boyle November 2010
fh = figure('position', [300 250 500 450], 'resize', 'off', 'windowstyle', 'modal');

bcolor = [0.8 0.8 0.8];

maps = {'jet', 'hot', 'gray', 'bone', 'ColorSpiral'};
mapValue = strmatch(colormap, maps);
if isempty(mapValue)
    mapValue =1;
end
uicontrol(fh, 'style', 'text', 'position', [125 400 250 30], 'string', 'Advanced Options', 'fontsize', 14, 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'text', 'position', [15 350 200 20], 'string', 'Stimulus Files Location:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
pathBox = uicontrol(fh, 'style', 'edit', 'position', [15 325 225 25], 'string', stimPath);
uicontrol(fh, 'style', 'pushbutton', 'position', [240 325 100 25], 'string', 'Browse...', 'callback', @browseFun);

uicontrol(fh, 'style', 'text', 'position', [15 285 75 20], 'string', 'colormap', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
mapMenu = uicontrol(fh, 'style', 'popupmenu', 'position', [15 265 120 20], 'string', maps, 'value', mapValue);
invertBox = uicontrol(fh, 'style', 'checkbox', 'position',[15 230 120 20], 'string', 'invert colormap', 'backgroundcolor', bcolor, 'value', invertColor);

bgh = uibuttongroup('parent', fh,  'units', 'pixels','position', [160 175 300 125], 'backgroundcolor', bcolor, 'title', 'colorbar range', 'selectionchangefcn', @selectionFun);

defColorButton = uicontrol(bgh, 'style', 'radiobutton', 'string', 'Default colorbars', 'position', [10 80 175 20], 'backgroundcolor', bcolor);

manColorButton = uicontrol(bgh, 'style', 'radiobutton', 'string', 'Manual colorbar range', 'position', [10 50 175 20], 'backgroundcolor', bcolor);

minBox = uicontrol(bgh, 'style', 'edit', 'position', [50 25 50 20]);
maxBox = uicontrol(bgh, 'style', 'edit', 'position', [150 25 50 20]);

uicontrol(bgh, 'style', 'text', 'string', 'to', 'position', [115 25 20 20], 'backgroundcolor', bcolor);

if length(colorRange) == 2
    set(bgh, 'selectedobject', manColorButton);   
    set(minBox, 'string', colorRange(1));
    set(maxBox, 'string', colorRange(2));
else
    set(bgh, 'selectedobject', defColorButton);
    set(minBox, 'enable', 'off');
    set(maxBox, 'enable', 'off');
end
%matchColorButton = uicontrol(bgh, 'style', 'radiobutton', 'string','Matchcolorbar ranges', 'position', [10 150 175 20], 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'pushbutton', 'position', [275 10 90 30], 'string', 'OK', 'callback', @okFun);

uicontrol(fh, 'style', 'pushbutton', 'position', [375 10 90 30], 'string', 'Cancel', 'callback', 'close(gcf)' );

if isempty(txtPath)
    chkValue = 0;
    enabled = 'off';
else
    chkValue = 1;
    enabled = 'on';
end
txtCheckBox = uicontrol(fh, 'style', 'checkbox', 'position',[15 125 200 20], 'string', 'send output to text file', 'backgroundcolor', bcolor, 'callback', @checkFun, 'value', chkValue);
uicontrol(fh, 'style', 'text', 'position', [35 100 200 20], 'string', 'output location:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
txtPathBox = uicontrol(fh, 'style', 'edit', 'position', [35 75 225 25], 'string', txtPath, 'enable', enabled);
txtPathButton = uicontrol(fh, 'style', 'pushbutton', 'position', [260 75 100 25], 'string', 'Browse...', 'callback', @txtBrowseFun, 'enable', enabled);

uiwait(fh);

    function browseFun(hobject, eventdata)
        thisPath = uigetdir(stimPath, 'Find Stimulus Folder');
        set(pathBox, 'string', [thisPath filesep]);
    end

    function txtBrowseFun(hObj, eventdata)
        thisPath = uigetdir(txtPath, 'Output text files location');
        set(txtPathBox, 'string', [thisPath filesep]);
    end

    function selectionFun(hObj, eventdata)
        if eventdata.NewValue == defColorButton
            set(minBox, 'enable', 'off');
            set(maxBox, 'enable', 'off');
        elseif eventdata.NewValue == manColorButton
            set(minBox, 'enable', 'on');
            set(maxBox, 'enable', 'on'); 
        end
    end

    function checkFun(hObj, eventdata)
        if(get(hObj, 'value'))
            set(txtPathBox, 'enable', 'on');
            set(txtPathButton, 'enable', 'on');
        else
            set(txtPathBox, 'enable', 'off');
            set(txtPathButton, 'enable', 'off');
        end
    end

    function okFun(hobject, eventdata)
       stimPath = get(pathBox, 'string');
       if ~isempty(stimPath) && ~exist(stimPath, 'dir')
           errordlg('designated folder does not exist, please re-enter');
           return
       elseif ~isempty(stimPath)
           if ~isequal(stimPath(end), filesep)
               stimPath = [stimPath filesep];
           end
       end
       
       colormap = maps{get(mapMenu, 'value')};
       invertColor = get(invertBox, 'value');
       selectedRadio = get(bgh, 'selectedobject');
       if selectedRadio == manColorButton
           min = str2num(get(minBox, 'string'));
           max = str2num(get(maxBox, 'string'));
           colorRange = [min max];
           if length(colorRange) ~= 2
               colorRange = [];
           end
       else
           colorRange = [];
       end
       if(get(txtCheckBox, 'value'))
           txtPath = get(txtPathBox, 'string');
           if ~exist(txtPath, 'dir')
               errordlg('designated folder does not exist, please re-enter');
               return
           end
       else
           txtPath = '';
       end
       close(fh);
    end
end