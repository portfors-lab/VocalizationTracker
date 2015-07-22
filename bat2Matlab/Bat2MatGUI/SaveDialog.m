
function [savePath format resolution] = SaveDialog(savePath, format, resolution)

if isempty(savePath)
    rootPath = pwd;
else
    rootPath = savePath;
end

formats = {'jpeg', 'pdf', 'eps','fig', 'tiff', 'png'};
resolutions = {'150', '300', '600'};

resolution = num2str(resolution);
formatValue = strmatch(format, formats);
resolutionValue = strmatch(resolution, resolutions);

fh = figure('position', [500 500 375 225], 'resize', 'off', 'windowstyle', 'modal');

bcolor = [0.8 0.8 0.8];
uicontrol(fh, 'style', 'text', 'position', [50 185 250 30], 'string', 'Save Parameters', 'fontsize', 14, 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'text', 'position', [15 150 100 20], 'string', 'Save Location', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);

pathBox = uicontrol(fh, 'style', 'edit', 'position', [15 125 225 25], 'string', savePath);

uicontrol(fh, 'style', 'text', 'position', [15 95 75 20], 'string', 'format', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'text', 'position', [180 95 75 20], 'string', 'resolution', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'pushbutton', 'position', [240 125 100 25], 'string', 'Browse...', 'callback', @browseFun);

saveFormat = uicontrol(fh, 'style', 'popupmenu', 'position', [15 75 120 20], 'string', formats, 'value', formatValue);

saveResolution = uicontrol(fh, 'style', 'popupmenu', 'position', [180 75 120 20], 'string', resolutions, 'value', resolutionValue);

uicontrol(fh, 'style', 'pushbutton', 'position', [175 10 90 30], 'string', 'OK', 'callback', @okFun);

uicontrol(fh, 'style', 'pushbutton', 'position', [275 10 90 30], 'string', 'Cancel', 'callback', 'close(gcf)' );

uiwait(fh);

    function browseFun(hObject, eventdata)
        savePath = uigetdir(rootPath, 'Save Figures Location');
        set(pathBox, 'string', [savePath filesep]);
    end

    function okFun(hObject, eventdata)
        formatValue = get(saveFormat, 'value');
        format = formats{formatValue};
        resolutionValue = get(saveResolution, 'value');
        resolution = resolutions{resolutionValue};
        resolution = str2double(resolution);
        savePath = get(pathBox, 'string');
        if ~exist(savePath, 'dir')
            errordlg('designated folder does not exist, please re-enter');
            return
        end
        if ~isequal(savePath(end), filesep)
            savePath = [savePath filesep];
        end
        close(fh);
    end
end