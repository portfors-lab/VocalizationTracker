function [figureHandle tableHandle tableData] = ExploreTestData(experimentData)
% ExploreTestData: displays a table of the contents of batlab experiment data

%Modified Amy Boyle Jan 2011
columnHeaders = {'Number of Traces','Test Type','Test Class','Comment'};
columnWidths = {75, 100, 125, 200};

try
    numTests = size(experimentData.test,2);
    tableData = cell(numTests,4);
    for rowNum = 1:numTests
        tableData{rowNum,1} = num2str(size(experimentData.test(rowNum).trace,2));
        tableData{rowNum,2} = experimentData.test(rowNum).testtype;
        tableData{rowNum,3} = experimentData.test(rowNum).full_testtype;
        tableData{rowNum,4} = experimentData.test(rowNum).comment;
    end
catch
    errordlg('file test data incomplete')
    return
end
figureHandle = figure('Name',[' Data Summary for File ' experimentData.pst_filename],'NumberTitle','off');

if getversion < 7.6 %R2008a
    figurePos = getpixelposition(figureHandle,0);
    figurePos(1:2) = 0;
    tableHandle = uitable(tableData, columnHeaders, 'Position', figurePos);
else
    tableHandle = uitable(figureHandle, 'Data', tableData, 'ColumnName', columnHeaders, 'ColumnWidth', columnWidths, 'units', 'normalized', 'position', [0 0 1 1]);% 'CellSelectionCallback', @cellSelectFun);
end

%set up double click callback for table
uitablepeer = findjobj(figureHandle, '-nomenu', 'class', 'uitablepeer');
set(uitablepeer,'MouseClickedCallback',@MouseClickHandler)

    function MouseClickHandler(handle,cbData)
    % handle ~ java object UITablePeer
    % cbData ~ callback data for the MouseClickedCallback event

        switch get(cbData,'ClickCount')
            case 1
    %             SingleClickHandler(handles,cbData);
            case 2
                DoubleClickHandler(handle,cbData);
            otherwise
                % unhandled for now
        end
    end

    function DoubleClickHandler(hObj, evtData)
        testNumber = get(hObj, 'SelectedRow')+1; %java, so indexes starting with 0
        ExploreTraceData(experimentData, testNumber);
    end
end