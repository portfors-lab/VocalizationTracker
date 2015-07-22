function [frameHandle, listboxHandle] = ExploreTestData(experimentData)

% %Create a preferences structure for the desired experimental data
% prefs = GeneratePreferences('mouse', '543', '', '');
% experimentData = LoadExperimentData(prefs);

numTests = size(experimentData.test,2);

tableData = cell(numTests+1,4);
tableData{1,1} = 'Test Number';
tableData{1,2} = 'Number of Traces';
tableData{1,3} = 'Test Type';
tableData{1,4} = 'Test Class';
tableData{1,5} = 'Comment';
for rowNum = 1:numTests
    tableData{rowNum+1,1} = experimentData.test(rowNum).testnum;
    tableData{rowNum+1,2} = size(experimentData.test(rowNum).trace,2);
    tableData{rowNum+1,3} = experimentData.test(rowNum).testtype;
    tableData{rowNum+1,4} = experimentData.test(rowNum).full_testtype;
    tableData{rowNum+1,5} = experimentData.test(rowNum).comment;
end

[frameHandle, listboxHandle] = gui_sheet('goo',tableData);
