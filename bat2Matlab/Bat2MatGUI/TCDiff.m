function TCDiff

%function TCD : an interface for finding the difference between two tuning
%curves.
%
% choose the folder in which your .pst and .raw files reside,
% this folder must have the SAME NAME as the .pst and .raw files
% you must also choose a range of test numbers over which to do the
% correlation. 
%
%The difference the absoulute value of the difference between two tuning
%curves, element-for-element.

%Amy Boyle Jan 2011

animalPath = [];
rootPath = pwd;
outPath = '';

tableIcon = imread('table.jpg');
%tableIcon = imresize(tableIcon, 0.8);

fh = figure('position', [300 250 425 450], 'resize', 'off', 'windowstyle', 'normal');
%fh = figure('position', [500 500 375 225], 'resize', 'off', 'windowstyle', 'modal');
bcolor = [0.8 0.8 0.8];

uicontrol(fh, 'style', 'text', 'position',  [25 400 350 30], 'string', 'Difference between two tuning curves', 'fontsize', 14, 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'text', 'position', [15 350 200 20], 'string', 'Animal folder', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);

pathBox = uicontrol(fh, 'style', 'edit', 'position', [15 325 225 25], 'string', animalPath);
uicontrol(fh, 'style', 'pushbutton', 'position', [240 325 100 25], 'string', 'Browse...', 'callback', @browseFun);

uicontrol(fh, 'style', 'pushbutton', 'position', [270 275 25 25],'CData', tableIcon, 'callback', @exploreFun);

uicontrol(fh, 'style', 'text', 'position', [15 275 100 20], 'string', 'Test numbers:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
minBox = uicontrol(fh, 'style', 'edit', 'position', [115 275 50 20]);
maxBox = uicontrol(fh, 'style', 'edit', 'position', [200 275 50 20]);

uicontrol(fh, 'style', 'text', 'string', '&', 'position', [170 275 20 20], 'backgroundcolor', bcolor);

% uicontrol(fh, 'style', 'text', 'position', [300 275 100 20], 'string', 'Trace:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
% traceBox = uicontrol(fh, 'style', 'edit', 'position', [350 275 50 20]);

uicontrol(fh, 'style', 'text', 'position', [15 225 100 20], 'string', 'Threshold:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
threshBox = uicontrol(fh, 'style', 'edit', 'position', [115 225 50 20], 'string', '0.2');

% txtCheckBox = uicontrol(fh, 'style', 'checkbox', 'position',[15 125 250 20], 'string', 'output data matrix to text file', 'backgroundcolor', bcolor, 'callback', @checkFun);
% uicontrol(fh, 'style', 'text', 'position', [35 100 200 20], 'string', 'output location:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
% txtPathBox = uicontrol(fh, 'style', 'edit', 'position', [35 75 225 25], 'string', outPath, 'enable', 'off');
% txtPathButton = uicontrol(fh, 'style', 'pushbutton', 'position', [260 75 100 25], 'string', 'Browse...', 'callback', @txtBrowseFun, 'enable', 'off');

uicontrol(fh, 'style', 'pushbutton', 'position', [200 10 90 30], 'string', 'OK', 'callback', @okFun);

uicontrol(fh, 'style', 'pushbutton', 'position', [300 10 90 30], 'string', 'Close', 'callback', 'close(gcf)');

prefs =[];
experiment_data = [];

    function browseFun(jObj,evtdata)
        animalPath = uigetdir(rootPath, 'folder containing data');
        if animalPath ~= 0
           set(pathBox, 'string', [animalPath filesep]);
           experiment_data = [];
        end
    end

    function exploreFun(hObj, eventData)
        %view file contents, get test numbers
        if isempty(animalPath)
            errordlg('! Must select animal file first, before able to select test numbers')
            return
        end
        if isempty(experiment_data)
            [prefs experiment_data] = GetExpData(animalPath);
            if isempty(prefs)
                return
            end
        end
        %ExploreTestData(experiment_data);
        testNums = SelectTestNumbers(experiment_data);
        if ~isempty(testNums)
            if length(testNums) ~= 2
                warndlg('you must select exactly 2 tests to find difference')
            else
                set(minBox, 'string', num2str(testNums(1)));
                set(maxBox, 'string', num2str(testNums(2)));
            end
        end
    end

    function okFun(oh,evtd)
        test1 = str2double(get(minBox, 'string'));
        test2 = str2double(get(maxBox, 'string'));
    
        animalPath = get(pathBox, 'string');
        if isempty(animalPath)
            disp('! You must enter a filepath for the data you wish to compare');
            return
        end            

        if isempty(experiment_data)
            [prefs experiment_data] = GetExpData(animalPath);
            if isempty(prefs)
                return
            end
        end
        prefs.spike_time_peak_threshold = str2double(get(threshBox, 'string'));

        if isnan(test1) || isnan(test2)
            disp('! Test range must inlcude numbers');
            return
        end
        if any([test1 test2] > length(experiment_data.test))
            disp('! Entered test range exceeds file contents')
            return
        end
        
        TuningCurveDifference(experiment_data, prefs, test1, test2);
    end
end