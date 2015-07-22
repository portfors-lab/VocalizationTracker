function VisualizationMachine(dirPath)
% function VisualizatonMachine
%
% An interface for graphing Batlab produced data
%
% make sure your javaclasspath is set for this to run
%
%This program will allow you to visually display batlab data.  Whole test
%analyses may be run to produce tuning curves or plot the PSTHs.
%Individual trace analysis allows several different plotting options,
%gathered on the same figure.

%close all
jpath=which('filetree.jar');
if all(cellfun(@isempty, strfind(javaclasspath, jpath))) %found filetree.jar in javaclasspath
    disp('adding filetree.jar to dynamic jpath')
    javaaddpath(jpath); %add it if not there
end

import java.awt.event.ActionEvent;
import filetree.*;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;
import javax.swing.JTextPane;
% import javax.swing.text.StyledDocument;
import javax.swing.JScrollPane;

%whether to use filesystem roots for the root of tree, if set to false it
%uses the user home
useRoots = false;

toolboxes = ver;
spbox = false;
statbox = false;
for toolbox = toolboxes;
    spbox = spbox || ~isempty(strfind(toolbox.Name, 'Signal Processing Toolbox'));
    statbox = statbox || ~isempty(strfind(toolbox.Name, 'Statistics Toolbox'));
end

if ~(spbox && statbox)
    errordlg('Missing toolbox: must have statistics and signal processing toolboxes installed');
    return;
end

err = 'red';
warn = 'orange';

if exist('dirPath', 'var')
    if exist(dirPath, 'dir')
        rootPath = dirPath;
    else
        errordlg(['Folder ' dirPath ' not found'])
        return
    end
else
    if useRoots
        rootPath
        h = [];
    else
        %os compatability
        if isunix || ismac
            rootPath = getenv('HOME');
        elseif ispc
            rootPath = getenv('USERPROFILE');
        else 
            warning('unrecognised/unsupported operating system');
            rootPath = pwd;
        end
        if isempty(rootPath)
            warning('home folder not found, using current directory for search path');
            rootPath = pwd;
        end
    end
end
%make sure the path ends with a slash
if ~isempty(rootPath)
if ~isequal(rootPath(end), filesep)
   rootPath = [rootPath filesep];
end
end

%default parameter values, load previous used and stored values

savePath = [pwd filesep];
format = 'tiff';
resolution = 150;
stimPath = [pwd filesep];
colormap = 'jet';
colorRange = [];
txtPath = [];
invertColor = 0;

prefPath = [rootPath 'machinePrefs.mat'];
if exist(prefPath, 'file')
    load(prefPath); %loads above declared variables
end
experiment_data = [];
prefs = [];
animalPath = [];

%==========================================================================

%this is the matlab portion of the figure, the java generated portion sits
%on top of this. the width of the javapanel will to match the width of the
%matlab figure
fh = figure;
set(fh, 'units', 'normalized','position', [0 0 1 0.3], 'menubar', 'none', 'name', 'Visualization Machine', 'numbertitle', 'off')
%set(fh, 'position', [150 150 900 100]);

SetMenuBar(fh, stimPath);

%main interface panel written in java
%create the java object
try
    if isempty(rootPath)
        javaPanel = DataSelectionPanel();
    else
        javaPanel = DataSelectionPanel(rootPath);
    end
catch e
    disp(['error:' e.message])
    disp('')
    disp('If you are using a version of matlab prior to R2010b, please call the GUI again, it should work this time')
    close(fh)
    return
end
%put the java object on the event dispatch thread
javaObjectEDT(javaPanel);
%display java object on the screen
[jObj container] = javacomponent(javaPanel, 'North', fh);
set(container, 'units', 'normalized', 'position', [0 0.1 1 0.8]); %position relative to matlab container
%set(container, 'position', [0 100 900 300]);

%popup menu for investigating what is in files before running graphs, it is
%necessary to do this in java as well, uicontextmenu does not work
menuItem = JMenuItem('File Contents...');
menuItem2 = JMenuItem('Preview spikes...');

set(menuItem, 'ActionPerformedCallback', @exploreFun);
set(menuItem2, 'ActionPerformedCallback', @previewFun);

jMenu = JPopupMenu;
jMenu.add(menuItem);
jMenu.add(menuItem2);

%==========================================================================
%matlab uicontrols

loadButton = uicontrol(fh, 'style', 'pushbutton', 'position', [15 30 150 50], 'string', 'Let''s Rock', 'callback', @loadFun);

panel = uipanel('parent', fh, 'units', 'pixels', 'position', [200 15 375 80], 'backgroundcolor', [0.8 0.8 0.8]);
uicontrol(panel, 'style', 'pushbutton', 'position', [10 15 140 25], 'string', 'Save Options...', 'callback', @saveOptions);
saveOption = uicontrol(panel, 'style', 'checkbox', 'position', [10 45 100 30], 'string', 'save figures', 'backgroundcolor', [0.8 0.8 0.8]);
uicontrol(panel, 'style', 'pushbutton', 'position', [200 15 140 25], 'string', 'Advanced Options...', 'callback', @advOptFun);

set(fh, 'CloseRequestFcn', @closeFun);

%need this info inorder to set doubdle clicking of java text fields-
%unfortunale result is tight coupling...
jvp = findjobj(javaPanel,'nomenu','class','javax.swing.JViewport'); 
modPanel = get(jvp, 'components');
% In order to eliminate a problem that occured with the use of the findjobj
% function, set the callback using the matlab handle, obtained from the
% handle function instead of setting it directly on the java object.
modPanel_mhandle = handle(modPanel(1), 'callbackProperties');
set(modPanel_mhandle, 'ComponentAddedCallback', @moduleAdded);

jTree = javaPanel.getTree();
jTree_mhandle = handle(jTree, 'callbackProperties');
set(jTree_mhandle, 'MousePressedCallback', {@mousePressedCallback,jMenu});

%status area

statusBox = JTextPane();
scrollPane = JScrollPane(statusBox);
statusText = ColorDocument();
statusBox.setStyledDocument(statusText);
[jhandle sbContainer] = javacomponent(scrollPane, [800 10 400 200], fh);
set(sbContainer, 'units', 'normalized', 'position', [0.60 0.05 0.38 0.9]);

    function saveOptions(hobject, eventdata)
        %allows user to select save parameters
        [savePath format resolution] = SaveDialog(savePath, format, resolution);
    end

    function advOptFun(hobject, eventdata)
        %allows user to select other parameters
        [stimPath colormap invertColor colorRange txtPath] = OptionDialog(stimPath, colormap, invertColor, colorRange, txtPath);
    end

    function loadFun(hobject, eventdata)
        %where the magic happens
%         set(hobject, 'enable', 'off');
        set(hobject, 'string', 'And Roll');
        drawnow;
        statusText.putLine('----------------Loading------------------');
        %function call here
        saveOn = get(saveOption, 'value');
        if saveOn
            RunVisualization(javaPanel, 'rootpath', rootPath, 'colorRange', colorRange, 'stimpath', stimPath, 'colormap', colormap, 'invertcolor', invertColor, 'saveon', {savePath format resolution});
        else
            RunVisualization(javaPanel, 'rootpath', rootPath, 'colorRange', colorRange, 'stimpath', stimPath, 'colormap', colormap, 'invertcolor', invertColor, 'txtpath', txtPath);
        end
         set(hobject, 'enable', 'on');
         set(hobject, 'string', 'Lets''s Rock');
        statusText.putLine('----------------Complete-----------------');
    end

    function mousePressedCallback(htree, eventData, jmenu)
        if eventData.isMetaDown
            x = eventData.getX;
            y = eventData.getY;
            jMenu.show(jTree, x, y);
            jMenu.repaint;
        end
    end

    function exploreFun(hObj, eventData)
        path = char(javaPanel.getCurrentSelection());
        if isempty(path)
            WriteStatus('Must highlight selection to view contents', err)
            return
        end
        path = [rootPath path];
        if isequal(path, animalPath)
            if isempty(experiment_data)
                [prefs experiment_data] = GetExpData(path);
            end
        else
            [prefs experiment_data] = GetExpData(path);
            animalPath = path;
        end
        if ~isempty(experiment_data)
            ExploreTestData(experiment_data);
        end
    end

    function previewFun(jObj, eventData)
        path = char(javaPanel.getCurrentSelection());
        if isempty(path)
            WriteStatus('Must highlight selection to view contents', err)
            return
        end
        path = [rootPath path];
        if isequal(path, animalPath)
            if isempty(experiment_data)
                [prefs experiment_data] = GetExpData(path);
            end
        else
            [prefs experiment_data] = GetExpData(path);
            animalPath = path;
        end
        if ~isempty(experiment_data)
            TraceChooser(prefs, experiment_data);
        end
    end

    function moduleAdded(hObj, eventData)
%         mods = findjobj(javaPanel,'nomenu','class', 'InputModule')
        mods = get(modPanel(1), 'components');
%         if ~isempty(mod)
%         for mod = mods
        for ind= 1:length(mods)
%             set(mod, 'MouseClickedCallback', @moduleClick)
%             [jO con] = javacomponent(javaPanel, 'North', fh);
            textFields = findjobj(mods(ind),'nomenu','class', 'javax.swing.JTextField');
            set(textFields, 'MouseClickedCallback', {@moduleClick, mods(ind)})
        end
    end
    
    function moduleClick(hObj, evtData, module)
        if get(evtData, 'clickCount') == 2
            name = get(module, 'name');
            paths = cellstr(char(javaPanel.getCheckedPaths()));
            indx = find(~cellfun('isempty', strfind(paths, [name filesep])));
            path = [rootPath paths{indx}];
            [prefs experiment_data] = GetExpData(path);
            switch hObj.name
                case 'test1'
                    if ~isempty(experiment_data)
                        nums = SelectTestNumbers(experiment_data);
                        if ~isempty(nums)
                            nums = num2str(nums');
                            javaPanel.setTests1(name, nums);
                        end
                    end
                    return
                case 'test2'
                    if ~isempty(experiment_data)
                        nums = SelectTestNumbers(experiment_data);
                        if ~isempty(nums)
                            nums = num2str(nums');
                            javaPanel.setTests2(name, nums);
                        end
                    end
                    return
                case 'trace'
                    tnums = str2num(javaPanel.getTests2(name));
                    switch char(javaPanel.getSelectionType(name));
                        case 'all'
                            return
                            %do nothing
                        case 'same'
                            if isempty(tnums)
                                WriteStatus('Select test numbers first', err);
                                return
                            end
                            traceCounts = [];
                            for t = tnums
                                traceCounts = [size(experiment_data.test(t).trace,2) traceCounts];
                            end
                            maxtrace = min(traceCounts);
                            liststr = {};
                            for tracenum = 1:maxtrace
                                liststr = [liststr {num2str(tracenum)}];
                            end
                            traceNums = listdlg('liststring', liststr, 'listsize', [150 150]);
                            if isempty(traceNums)
                                return
                            end
                            traceNums = num2str(traceNums);
                        case 'different'
                            if isempty(tnums)
                                WriteStatus('Select test numbers first', err);
                                return
                            end
                            traceNums = [];
                            for t = tnums
                                nums = num2str(SelectTraceNumbers(experiment_data, t)');
                                traceNums = [traceNums '; ' nums];
                                if isempty(nums)
                                    return
                                end
                            end
                            traceNums = traceNums(2:end); %remove leading ;
                        otherwise
                            return
                    end
                    javaPanel.setTraces(name, traceNums)
            end
        end
    end

    function closeFun(jObj, eventData)
        try
            save(prefPath, 'savePath', 'stimPath', 'format', 'resolution', 'colormap', 'invertColor', 'colorRange', 'txtPath');
        catch e
            disp('Warning: did not save GUI user preferences');
            disp(e.message);
        end
        delete(fh);        
    end
        
end