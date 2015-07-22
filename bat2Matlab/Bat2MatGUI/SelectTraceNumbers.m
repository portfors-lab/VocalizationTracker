function [traceNums] = SelectTraceNumbers(experimentData, testNum)

% %Create a preferences structure for the desired experimental data
% prefs = GeneratePreferences('mouse', '543', '', '');
% experimentData = LoadExperimentData(prefs);

% Last modified Amy Boyle 5/2/11, added multistimulus support

numTraces = size(experimentData.test(testNum).trace,2);
testtype = experimentData.test(testNum).testtype;

columnHeaders{1,1} = 'Record Duration';
columnHeaders{1,2} = 'Number of Sweeps';
columnHeaders{1,3} = 'Attenuation';
columnHeaders{1,4} = 'Stimulus Duration';
columnHeaders{1,5} = 'Delay';
% if strcmp(testtype,'tone')
%     columnHeaders{1,6} = 'Frequency';
% elseif strcmp(testtype,'twotone')
%     columnHeaders{1,6} = 'Frequency 1';
%     columnHeaders{1,7} = 'Frequency 2';
if strcmp(testtype,'vocalization')
    columnHeaders{1,6} = 'File Name';
else
    columnHeaders{1,6} = 'Frequency';
end

tableData = cell(numTraces,size(columnHeaders,2));

for rowNum = 1:numTraces
    tableData{rowNum,1} = num2str(experimentData.test(testNum).trace(rowNum).record_duration);
    tableData{rowNum,2} = num2str(experimentData.test(testNum).trace(rowNum).num_samples);
    if ~experimentData.test(testNum).trace(rowNum).is_control
        if length(experimentData.test(testNum).trace(rowNum).stimulus) > 1
            stimCount = length(experimentData.test(testNum).trace(rowNum).stimulus);
            attenuationArr = zeros(1,stimCount);
            durationArr = zeros(1,stimCount);
            delayArr = zeros(1,stimCount);
            frequencyArr = zeros(1,stimCount);
            vocalArr = [];
            for stimNum = 1:stimCount
                attenuationArr(stimNum) = experimentData.test(testNum).trace(rowNum).stimulus(stimNum).attenuation;
                durationArr(stimNum) = experimentData.test(testNum).trace(rowNum).stimulus(stimNum).duration;
                delayArr(stimNum) = experimentData.test(testNum).trace(rowNum).stimulus(stimNum).delay;
                frequencyArr(stimNum) = experimentData.test(testNum).trace(rowNum).stimulus(stimNum).frequency;
                vocalArr = [vocalArr ', ' experimentData.test(testNum).trace(rowNum).stimulus(stimNum).vocal_call_file];
            end
            if (sum(attenuationArr/attenuationArr(1))) == stimCount
                attenuationArr = attenuationArr(1);
            end
            if sum(durationArr/durationArr(1)) == stimCount
                durationArr = durationArr(1);
            end
            if sum(delayArr/durationArr(1)) == stimCount
                delayArr = delayArr(1);
            end
            if sum(frequencyArr/frequencyArr(1)) == stimCount
                frequencyArr = frequencyArr(1);
            end
            tableData{rowNum,3} = num2str(attenuationArr);
            tableData{rowNum,4} = num2str(durationArr);
            tableData{rowNum,5} = num2str(delayArr);
            if strcmp(testtype,'vocalization')
                tableData{rowNum,6} = vocalArr;
            else
                tableData{rowNum,6} = num2str(frequencyArr);
            end
        else
            tableData{rowNum,3} = num2str(experimentData.test(testNum).trace(rowNum).stimulus.attenuation);
            tableData{rowNum,4} = num2str(experimentData.test(testNum).trace(rowNum).stimulus.duration);
            tableData{rowNum,5} = num2str(experimentData.test(testNum).trace(rowNum).stimulus.delay);
            if strcmp(testtype,'tone')
                tableData{rowNum,6} = num2str(experimentData.test(testNum).trace(rowNum).stimulus.frequency);
            elseif strcmp(testtype,'twotone')
                tableData{rowNum,6} = num2str(experimentData.test(testNum).trace(rowNum).stimulus.frequency);
                tableData{rowNum,7} = num2str(experimentData.test(testNum).trace(rowNum).stimulus(2).frequency);
            elseif strcmp(testtype,'vocalization')
                tableData{rowNum,6} = experimentData.test(testNum).trace(rowNum).stimulus.vocal_call_file;
            end
        end
    else
        tableData{rowNum,3} = 'N/A';
        tableData{rowNum,4} ='N/A';
        tableData{rowNum,5} ='N/A';
        if strcmp(testtype,'tone')
            tableData{rowNum,6} = 'N/A';
        elseif strcmp(testtype,'twotone')
            tableData{rowNum,6} = 'N/A';
            tableData{rowNum,7} = 'N/A';
        elseif strcmp(testtype,'vocalization')
            tableData{rowNum,6} = 'N/A';
        end
    end
end

figureHandle = figure('Name',[' Data Summary for Test ' int2str(testNum)],'NumberTitle','off');

figurePos = getpixelposition(figureHandle,1);
figurePos(1:2) = 0;

tableHandle = uitable(figureHandle, 'Data',tableData, 'ColumnName',columnHeaders, 'units', 'normalized', 'Position', [0 0 1 1]);
uicontrol(figureHandle, 'style', 'pushbutton', 'position', [200 10 75 25], 'string', 'OK', 'callback', 'close(gcf)');

uicontrol(figureHandle, 'style', 'pushbutton', 'position', [300 10 75 25], 'string', 'Cancel', 'callback', @cancelFun);
traceNums  =[];

%set up double click callback for table
set(tableHandle, 'cellselectioncallback', @cellSelectFun);

uiwait(figureHandle);

    function cellSelectFun(hObj, eventData)
     traceNums = eventData.Indices(:,1);
    end

    function cancelFun(hObj, eventData)
        traceNums = [];
        close(figureHandle)
    end
end