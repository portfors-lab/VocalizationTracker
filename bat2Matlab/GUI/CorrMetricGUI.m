function CorrMetricGUI
% function CorrMetricGUI  :  An interface for correlating Batlab data
% 
% choose the folder in which your .pst and .raw files reside,
% this folder must have the SAME NAME as the .pst and .raw files
% you must also choose a range of test numbers over which to do the
% correlation. Select the table icon to see file contents, then select the
% desired columns using the shift and control keys.  To see the raw data
% sweeps select the trace sweep icon in the top right corner of the window
%
% This function will produce 4 graphs, which are:
% 1: figure that plots the selected trace from each test
% 2: figure that plots the first sweep of selected trace from each test
% 3: color plot of the correlation between the tests
%
%Data output to a text file is titled the animal number.
%In the file, at the top of the columns are the names of the stimulus 
%used for that test, these same names will then also correspond down the 
%rows starting top to bottom. i.e. column 1 name == row 1 name.
%

%Amy Boyle Jan 2011, adapted from testCorrMetric.m
%close all
clear variables

%%%%%%%%%%%%%%%%%%%%%%%
%GUI

animalPath = [];
rootPath = pwd;
outPath = '';

tableIcon = imread('table.jpg');
%tableIcon = imresize(tableIcon, 0.8);
traceIcon = imread('trace.jpg');
%traceIcon = imresize(traceIcon, 0.05);

fh = figure('position', [300 250 425 450], 'resize', 'off', 'windowstyle', 'normal');
%fh = figure('position', [500 500 375 225], 'resize', 'off', 'windowstyle', 'modal');
bcolor = [0.8 0.8 0.8];

uicontrol(fh, 'style', 'text', 'position',  [100 400 250 30], 'string', 'Correlation Matrix', 'fontsize', 14, 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'text', 'position', [15 350 200 20], 'string', 'Animal folder', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);

pathBox = uicontrol(fh, 'style', 'edit', 'position', [15 325 225 25], 'string', animalPath);
uicontrol(fh, 'style', 'pushbutton', 'position', [240 325 100 25], 'string', 'Browse...', 'callback', @browseFun);

uicontrol(fh, 'style', 'pushbutton', 'position', [115 275 25 25],'CData', tableIcon, 'tooltipstring', 'select test numbers', 'callback', @exploreFun);
uicontrol(fh, 'style', 'pushbutton', 'position', [350 325 25 25],'CData', traceIcon, 'tooltipstring', 'view spike traces',  'callback', @previewFun);

uicontrol(fh, 'style', 'text', 'position', [15 275 100 20], 'string', 'Test numbers:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
testNumBox =  uicontrol(fh, 'style', 'edit', 'position', [150 275 150 20], 'callback', @resetNums);

%right now, this code just uses first trace, however, if desired, choosing
% the trace number can be enabled by uncommenting the following
uicontrol(fh, 'style', 'text', 'position', [310 275 100 20], 'string', 'Trace:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
traceBox = uicontrol(fh, 'style', 'edit', 'position', [350 275 50 20], 'string', '1');

bgh = uibuttongroup('parent', fh,  'units', 'pixels','position', [15 175 400 75], 'backgroundcolor', bcolor, 'selectionchangefcn', @selectionFun);

universalRadio = uicontrol(bgh, 'style', 'radiobutton', 'string', 'Universal Threshold', 'position', [10 40 175 20], 'backgroundcolor', bcolor);
independentRadio = uicontrol(bgh, 'style', 'radiobutton', 'string', 'Independent Thresholds', 'position', [175 40 200 20], 'backgroundcolor', bcolor);
threshBox = uicontrol(bgh, 'style', 'edit', 'position', [25 10 50 20], 'string', '0.2');
indThreshButton = uicontrol(bgh, 'style', 'pushbutton', 'position', [175 10 200 25], 'string', 'independent thresholds', 'enable', 'off', 'callback', @threshFun);

txtCheckBox = uicontrol(fh, 'style', 'checkbox', 'position',[15 125 250 20], 'string', 'output data matrix to text file', 'backgroundcolor', bcolor, 'callback', @checkFun);
uicontrol(fh, 'style', 'text', 'position', [35 100 200 20], 'string', 'output location:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
txtPathBox = uicontrol(fh, 'style', 'edit', 'position', [35 75 225 25], 'string', outPath, 'enable', 'off');
txtPathButton = uicontrol(fh, 'style', 'pushbutton', 'position', [260 75 100 25], 'string', 'Browse...', 'callback', @txtBrowseFun, 'enable', 'off');

okButton = uicontrol(fh, 'style', 'pushbutton', 'position', [200 10 90 30], 'string', 'OK', 'callback', @okFun);

uicontrol(fh, 'style', 'pushbutton', 'position', [300 10 90 30], 'string', 'Close', 'callback', 'close(gcf)');

prefs = [];
experiment_data = [];
thresholds = [];

    function browseFun(jObj,evtdata)
        animalPath = uigetdir(rootPath, 'folder containing data');
        if animalPath ~= 0
           set(pathBox, 'string', [animalPath filesep]);
           experiment_data = [];
           thresholds = [];
           set(testNumBox, 'string', '');
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
        set(testNumBox, 'string', num2str(testNums'));
        thresholds = [];
    end

    function previewFun(hObj, eventData)
        %view individual trace sweep
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
        traceChooser(prefs, experiment_data);
    end

    function txtBrowseFun(hObj, eventdata)
        thisPath = uigetdir(rootPath, 'Output text files location');
        set(txtPathBox, 'string', [thisPath filesep]);
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

    function selectionFun(hObj, eventdata)
        if eventdata.NewValue == universalRadio
           set(indThreshButton, 'enable', 'off');
           set(threshBox, 'enable', 'on');
        elseif eventdata.NewValue == independentRadio
           set(threshBox, 'enable', 'off');
           set(indThreshButton, 'enable', 'on');
        end
    end

    function threshFun(hObj, eventData)
        testNums = str2num(get(testNumBox, 'string'));
        if isempty(testNums)
            errordlg('! Must select tests to correlate');
            return
        end
        thresholds = ThresholdDialog(testNums, thresholds);
    end

    function resetNums(hObj, evtData)
        thresholds = [];
    end

    function okFun(hObj,eventData)
        %set(okButton, 'enable', 'off'); --when not commented out, this
        %executes out of order and doesn't disable until right before
        %algorigthm finishes.
        
        calculate_spike_rates = true;
        temporalSimilarity = true;        
        
        %Get input values and error check------------------
        test_nums = str2num(get(testNumBox, 'string'));
        if isempty(test_nums)
            errordlg('! Must select tests to correleate');
            return
        end
        outputMat = get(txtCheckBox, 'value');
        if outputMat
            outPath = get(txtPathBox, 'string');
            if isempty(outPath)
                errordlg('! You must specify a path to save to')
                return
            end
            if ~isequal(outPath(end), filesep)
                outPath = [outPath filesep];
            end
        end
        
        animalPath = get(pathBox, 'string');
        if isempty(animalPath)
            errordlg('! You must enter a filepath for the data you wish to correlate');
            return
        end

        if exist('traceBox')
            trace_num = str2double(get(traceBox, 'string'));
            if isnan(trace_num)
                trace_num = 1;
                disp('! trace number not recognized, using trace 1');
            end
        else
            trace_num = 1;
        end

        %load the experimental data
        if isempty(experiment_data)
            [prefs experiment_data] = GetExpData(animalPath);
            if isempty(prefs)
                return
            end
        end        
        
        % done error checking------------------------------          
               
        selectedRadio = get(bgh, 'selectedobject');
        if selectedRadio == universalRadio
           prefs.spike_time_peak_threshold = str2double(get(threshBox, 'string'));
           thresholds = [];
        else
           if isempty(thresholds) 
               thresholds = ones(length(test_nums))*str2double(get(threshBox, 'string'));
           end
        end
        
        %get the stimulus names for each test to use for plot
         count = 1;
         for thisTest = test_nums
           if isempty(experiment_data.test(thisTest).trace(trace_num).stimulus)
               vocalName = experiment_data.test(thisTest).testtype;
           else
               vocalName = experiment_data.test(thisTest).trace(trace_num).stimulus(1).vocal_call_file;
               
               if isempty(vocalName)
                   vocalName = experiment_data.test(thisTest).testtype;
               else
                   vocalName = vocalName(1:end-6);
               end
               vocalName = [vocalName '-' num2str(thisTest)];
           end
           stimulusNames{count} = vocalName;
           count = count + 1;
         end        
        num_tests = length(test_nums);
        if calculate_spike_rates
            display('Calculating Spike Rates');
            spikeTimes = CalculateSpikeTimes(experiment_data, prefs, test_nums, trace_num, [], thresholds);
            % FilterSpikeTrain uses a Gaussian Kernel for smoothing
            spike_rates = FilterSpikeTrain(experiment_data, ...
                                             prefs, ....
                                             test_nums, ...
                                             trace_num, ...
                                             spikeTimes);
           %save ('testCorrMetrixSpikeRates', 'spike_rates');
        else
            display('Loading spike Rates off Disk');
            load testCorrMetrixSpikeRates
        end
        
        display('Generating Distance Matrix');
        tempSim_correlation_matrix = zeros(num_tests,num_tests);
        num_plot_rows = ceil(sqrt(num_tests));
        num_plot_cols = num_plot_rows;
        test_num_idxs = 1:length(test_nums);
        for test_num1_idx = test_num_idxs
            test_num1 = test_nums(test_num1_idx);
            for test_num2_idx = 1:test_num1_idx;
                test_num2 = test_nums(test_num2_idx);
%                 if test_num1 == test_num2
%                     consistancyCorrelation1 = CalculateSpikeCorr(experiment_data, ...
%                                                        prefs, ...
%                                                        test_num1, ...
%                                                        trace_num, ...
%                                                        test_num1, ...
%                                                        trace_num, ...
%                                                        spike_rates...
%                                                        ); 
%                                                    tempSim_correlation_matrix(test_num1_idx, test_num2_idx) = consistancyCorrelation1;
% 
%                 else
                correlation = CalculateSpikeCorr2(experiment_data, ...
                                                   prefs, ...
                                                   test_num1, ...
                                                   trace_num, ...
                                                   test_num2, ...
                                                   trace_num, ...
                                                   spike_rates);
                   consistancyCorrelation1 = CalculateSpikeCorr(experiment_data, ...
                                                       prefs, ...
                                                       test_num1, ...
                                                       trace_num, ...
                                                       test_num1, ...
                                                       trace_num, ...
                                                       spike_rates...
                                                       );
                   consistancyCorrelation2 = CalculateSpikeCorr(experiment_data, ...
                                                       prefs, ...
                                                       test_num2, ...
                                                       trace_num, ...
                                                       test_num2, ...
                                                       trace_num, ...
                                                       spike_rates ...
                                                       );
                   meanDiff = CalculateRateDiff(experiment_data, ...
                                                       prefs, ...
                                                       test_num1, ...
                                                       trace_num, ...
                                                       test_num2, ...
                                                       trace_num, ...
                                                       spike_rates);
                   selectivity_index_matrix(test_num1_idx, test_num2_idx) = meanDiff;
                   selectivity_index_matrix(test_num2_idx, test_num1_idx) = meanDiff;
                   correlation_matrix(test_num1_idx, test_num2_idx) = correlation;
                   correlation_matrix(test_num2_idx, test_num1_idx) = correlation;
                   c1(test_num1_idx, test_num2_idx) = consistancyCorrelation1;
                   c1(test_num2_idx, test_num1_idx) = consistancyCorrelation1;
                   c2(test_num1_idx, test_num2_idx) = consistancyCorrelation2;
                   c2(test_num2_idx, test_num1_idx) = consistancyCorrelation2;
                   tempSim_correlation_matrix(test_num1_idx, test_num2_idx) = correlation/(consistancyCorrelation1 + consistancyCorrelation2);
                   tempSim_correlation_matrix(test_num2_idx, test_num1_idx) = correlation/(consistancyCorrelation1 + consistancyCorrelation2);
%                  end
            end
        end
        
        %Just plot the first trace of each test, to get a better sense for
        %individual correlations
        figure;
        test_num1_idx = 0;
        for test_num1 = test_nums
            test_num1_idx = test_num1_idx + 1;
            subplot(num_plot_rows,num_plot_cols,test_num1_idx);
            plot(spike_rates{test_num1,trace_num}(1,:));
            %title(['Test ' int2str(test_num1)]);
            title(stimulusNames{test_num1_idx});
        end

        figure;
        test_num1_idx = 0;
        for test_num1 = test_nums
            test_num1_idx = test_num1_idx + 1;
            subplot(num_plot_rows,num_plot_cols,test_num1_idx);
            plot(spike_rates{test_num1,trace_num}');
            %title(['Test ' int2str(test_num1)]);
            title(stimulusNames{test_num1_idx});
        end     
        
        figure;
        imagesc(test_num_idxs,test_num_idxs,selectivity_index_matrix);
        title(['Selectivity Index for ' prefs.cell_id4_plot]);
        set(gca,'YTick', test_num_idxs, 'YTickLabel', stimulusNames, 'XTick', test_num_idxs, 'XTickLabel', stimulusNames); %stimulus name labels
        xticklabel_rotate; %rotate xaxis labels
        colorbar
        
%         c1(find(c1 < 0)) = 0;
%         figure;
%         imagesc(test_num_idxs,test_num_idxs,c1);
%         title(['Consistency Correlation1 (Rx) for ' prefs.cell_id4_plot]);
%         set(gca,'YTick', test_num_idxs, 'YTickLabel', stimulusNames, 'XTick', test_num_idxs, 'XTickLabel', stimulusNames); %stimulus name labels
%         xticklabel_rotate; %rotate xaxis labels
%         colorbar
%         
%         c2(find(c2 < 0)) = 0;
%         figure
%         imagesc(test_num_idxs,test_num_idxs,c2);
%         title(['Consistency Correlation2 (Rx) for ' prefs.cell_id4_plot]);
%         set(gca,'YTick', test_num_idxs, 'YTickLabel', stimulusNames, 'XTick', test_num_idxs, 'XTickLabel', stimulusNames); %stimulus name labels
%         xticklabel_rotate; %rotate xaxis labels
%         colorbar
        
        correlation_matrix(find(correlation_matrix < 0)) = 0;
        figure;
        imagesc(test_num_idxs,test_num_idxs,correlation_matrix);
        title(['Mean Cross Correlation(Rxy) for ' prefs.cell_id4_plot]);
        set(gca,'YTick', test_num_idxs, 'YTickLabel', stimulusNames, 'XTick', test_num_idxs, 'XTickLabel', stimulusNames); %stimulus name labels
        xticklabel_rotate; %rotate xaxis labels
        colorbar
        
        tempSim_correlation_matrix(find(tempSim_correlation_matrix < 0)) = 0;
        figure;
        imagesc(test_num_idxs,test_num_idxs,tempSim_correlation_matrix);
        title(['Temporal Similarity Mean Correlation(Sxy) for ' prefs.cell_id4_plot]);
        set(gca,'YTick', test_num_idxs, 'YTickLabel', stimulusNames, 'XTick', test_num_idxs, 'XTickLabel', stimulusNames); %stimulus name labels
        xticklabel_rotate; %rotate xaxis labels
        colorbar
        %caxis([0.4 0.7])
        
        %confusion matrix, may use for any correlation measure
        confusionMatrix = MatrixInformationTransfer(tempSim_correlation_matrix);
        
        %if desired, output the correlation matrix to text file
        if outputMat
            fname1 = [outPath prefs.cell_id '_crossCorrelationData.txt'];
            corrFid1 = fopen(fname1, 'w');
            fname2 = [outPath prefs.cell_id '_tempSimCorrelationData.txt'];
            corrFid2 = fopen(fname2, 'w');
            fname3 = [outPath prefs.cell_id '_selectIndexData.txt'];
            corrFid3 = fopen(fname3, 'w');
            for name = stimulusNames %print stimulus name headers
                fprintf(corrFid1, '%s,', name{1});
                fprintf(corrFid2, '%s,', name{1});
                fprintf(corrFid3, '%s,', name{1});
            end
            fprintf(corrFid1, '\n');
            fprintf(corrFid2, '\n');
            fprintf(corrFid3, '\n');
            fclose(corrFid1);
            fclose(corrFid2);
            fclose(corrFid3);
            dlmwrite(fname1, correlation_matrix, '-append');
            dlmwrite(fname2, tempSim_correlation_matrix, '-append');
            dlmwrite(fname3, selectivity_index_matrix, '-append');
        end

% TODO: This is the cluster analysis that was used in previous version of the Correlation Metric code,
% the next step here is to implement this for the temporal similarity
% correlation
%
%         for i = 1:num_tests
%             test_num = test_nums(i);
%             this_distance = correlation_matrix(i,i);
%             this_stdev = stdev_correlation_matrix(i,i);
%             matches = [];
%             for j = 1:num_tests
%                 corr_distance = correlation_matrix(i,j);
%                 corr_stdev = stdev_correlation_matrix(i,j);
%                 if abs(this_distance - corr_distance) <= 0.3*(this_stdev + corr_stdev)
%                     matches = [matches test_nums(j)];
%                 end
%             end
%             display(['Cluster ' int2str(test_num) ':  ' num2str(matches)]);
%         end
%         
%         %%%This performs a cluster analysis of the correlation data.
%         correlation_matrix = 1 - correlation_matrix;
%         for i = 1:length(correlation_matrix)
%             correlation_matrix(i,i) = 0;
%         end
%         num_clusters = 5;
%         vector_correlation = squareform(correlation_matrix);
%         Z = linkage(vector_correlation);
%         T = cluster(Z,'maxclust',num_clusters);
%         % T = [test_nums' T];
%         for i = 1:num_clusters
%             cluster_traces = test_nums(find(T == i));
%             display(['Cluster ' int2str(i) ':  ' num2str(cluster_traces)]);
%         end
         set(okButton, 'enable', 'on')
     end     
end