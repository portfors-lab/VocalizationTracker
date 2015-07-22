function TextFileOutput(experiment_data, prefs, testNums, traceNumArray, flags, txtPath)

%     outputFolder = [prefs.bat2matlab_directory prefs.cell_id 'txtdata'];
    outputFolder = [txtPath prefs.cell_id 'txtdata'];
    if ~exist(outputFolder, 'file')
        mkdir(outputFolder);
    end
%need to combine cell array to get all trace nums for rawData
%!!! This will cause an index exceeds matrix dimensions error if
    if iscell(traceNumArray)
        totalTraces = cell2mat(traceNumArray');
        totalTraces = unique(totalTraces);
    else
        totalTraces = traceNumArray;
    end
    rawData = ExtractRawData(experiment_data, prefs, testNums, totalTraces);
    
    plot_signal             = flags(1);
    plot_smoothed_signal    = flags(2);
    plot_detected_spikes    = flags(3);
    plot_histogram          = flags(4);
    plot_smoothed_histogram = flags(5);
    plot_stimulus           = flags(6);
    plot_spike_rates        = flags(7);
    plot_spectrogram        = flags(8);
    plot_peaks              = flags(9);
    plot_dB_SPL             = flags(10);

    
%     if ~isequal(unique(testNums), testNums)
%         [testNums, ~, n] = unique(testNums);
%         temp = cell(size(testNums));
%         for ind = 1:length(traceNumArray);
%             temp{n(ind)} = [temp{n(ind)} traceNumArray{ind}];
%         end
%         traceNumArray = temp;
%     end
    try
        if plot_detected_spikes || plot_spike_rates
                    [spike_times spike_idxs smoothed_data peak_thresholds] ...
                        = CalculateSpikeTimes(experiment_data, ...
                                                        prefs, ...
                                                        testNums, ...
                                                        totalTraces, ...
                                                        rawData);
        end

        if plot_histogram || plot_smoothed_histogram
            [hists sweepHists binWidth binCenters] = ...
                GenerateHistograms(experiment_data, prefs, testNums, totalTraces);
        end

        if plot_stimulus || plot_dB_SPL
            [individual_stimulus_signals summed_stimulus_signals stimulus_sampling_frequencies] = ...
                GenerateStimulus(experiment_data, ...
                                           prefs, ...
                                           testNums, ...
                                           totalTraces);
        end

%         if plot_dB_SPL
%             [stimulus_spectrograms stimulus_spectrogram_time_idxs stimulus_spectrogram_freq_idxs stimulus_dB_SPLs] =  ...
%                 GenerateSpectrograms(experiment_data, ...
%                                                prefs, ...
%                                                testNums, ...
%                                                totalTraces, ...
%                                                summed_stimulus_signals);
%         end
        if plot_spectrogram || plot_peaks
            disp('! text output for spectrogram and peaks not supported');
        end
    catch e
        disp(['! ' e.message ]);     
    end
    if isempty(traceNumArray)
        allTraces = true;
    else
        allTraces = false;
    end
    traceListIndex = 1;
    for test = testNums
        if plot_signal
            signalFile = [outputFolder filesep 'test_' num2str(test) 'MER.txt'];
            signalFid = fopen(signalFile, 'w');
        end
        if plot_smoothed_signal
            smthsignalFile = [outputFolder filesep 'test_' num2str(test) 'smthMER.txt'];
            smthsignalFid = fopen(smthsignalFile, 'w');            
        end
        if plot_detected_spikes
            spikesFile = [outputFolder filesep 'test_' num2str(test) 'spikes.txt'];
            spikesFid = fopen(spikesFile, 'w');
        end
        if plot_histogram
            histsFile = [outputFolder filesep 'test_' num2str(test) 'hists.txt'];
            histsFid = fopen(histsFile, 'w');
        end
        if plot_smoothed_histogram
            smthhistsFile = [outputFolder filesep 'test_' num2str(test) 'smthhists.txt'];
            smthhistsFid = fopen(smthhistsFile, 'w');
        end
        if plot_stimulus
            stimFile = [outputFolder filesep 'test_' num2str(test) 'stim.txt'];
            stimFid = fopen(stimFile, 'w');
        end
        if plot_spike_rates
            spikeRateFile = [outputFolder filesep 'test_' num2str(test) 'spikerate.txt'];
            spikeRateFid = fopen(spikeRateFile, 'w'); 
        end
        if plot_dB_SPL
            dBFile = [outputFolder filesep 'test_' num2str(test) 'dBSPL.txt'];
            dBFid = fopen(dBFile, 'w');
        end
        if allTraces
            traceNums = 1:length(experiment_data.test(test).trace); 
        elseif iscell(traceNumArray)
            traceNums = traceNumArray{traceListIndex};
        else
            traceNums = traceNumArray;
        end
        index = 1;        
        merMat = [];
        smthmerMat =[];
        for trace = traceNums
            if plot_signal
                merMat = [merMat rawData{test,trace}.'];
                numSweeps = size(rawData{test,trace},1);
                for sweep = 1:numSweeps
                    fprintf(signalFid, 'trace%dr%d\t', trace, sweep);
                end
            end
            if plot_smoothed_signal
                smthmerMat = [smthmerMat smoothed_data{test,trace}.'];
                numSweeps = size(smoothed_data{test,trace},1);
                for sweep = 1:numSweeps
                    fprintf(smthsignalFid, 'trace%dr%d\t', trace, sweep);
                end
            end
            if plot_detected_spikes                    
                fprintf(spikesFid, 'trace%d\t', trace);
                detSpikes = spike_times{test,trace};
                detSpikes = cell2mat(detSpikes);
                spikeMat(1:length(detSpikes),index) = detSpikes;
            end
            if plot_histogram
               bc = binCenters{test,trace};
               histogram = hists{test,trace};
               %Normalize the histogram to average spikes per bin
               histogram = histogram / size(sweepHists{test,trace},1);
               histMat(:,index*2-1) = bc';
               histMat(:,index*2) = histogram';
               fprintf(histsFid, 'trace%dbinCen\ttrace%dhist\t', trace, trace);
            end
            if plot_smoothed_histogram
               bc = binCenters{test,trace};
               histogram = hists{test,trace};
               %Normalize the histogram to average spikes per bin
               histogram = histogram / size(sweepHists{test,trace},1);
               y = SmoothSeries(bc,histogram,bc,2);
               smthhistMat(:,index*2-1) = bc';
               smthhistMat(:,index*2) = y';
               fprintf(smthhistsFid, 'trace%dbinCen\ttrace%dhist\t',trace, trace);
            end
            if plot_stimulus
                [idx stim] = processStimulus(test,trace);
                stimMat(:, index*2-1) = idx';
                stimMat(:, index*2) = stim';
                fprintf(stimFid, 'trace%didx\ttrace%dstim\t', trace, trace);
            end
            if plot_spike_rates
                [spike_rates, spike_rate_averages, spike_rate_sampling_frequencies] ...
                    = GenerateSpikeRates(experiment_data, ...
                                           prefs, ...
                                           test, ...
                                           trace, ...
                                           spike_idxs);
                signal_duration = experiment_data.test(test).trace(trace).record_duration; %In milliseconds
                spike_rate = spike_rate_averages{test,trace};  
                spike_rate = spike_rate / 1000; %In spikes per millisecond
                spike_rate_idx = linspace(0,signal_duration,length(spike_rate));
                spikeRateMat(:,index*2-1) = spike_rate_idx';
                spikeRateMat(:, index*2) = spike_rate';
                fprintf(spikeRateFid, 'trace%dindx\ttrace%drate\t', trace, trace);
            end
            if plot_dB_SPL
                            [~, ~, ~, stimulus_dB_SPLs] =  ...
                GenerateSpectrograms(experiment_data, ...
                                               prefs, ...
                                               test, ...
                                               trace, ...
                                               summed_stimulus_signals);
                dBMat(:, index) = stimulus_dB_SPLs{test, trace}';
                fprintf(dBFid, 'trace%d\t', trace);
            end
            index = index+1;
            
        end
        traceListIndex = traceListIndex+1;  
        if plot_signal
            fprintf(signalFid, '\n');
            fclose(signalFid);
            dlmwrite(signalFile, merMat, '-append', 'delimiter', '\t');
            clear merMat;
        end
        if plot_smoothed_signal
            fprintf(smthsignalFid, '\n');
            fclose(smthsignalFid);
            dlmwrite(smthsignalFile, smthmerMat, '-append', 'delimiter', '\t');
            clear smthmerMat;
        end
        if plot_detected_spikes
            fprintf(spikesFid, '\n');
            fclose(spikesFid);
            dlmwrite(spikesFile, spikeMat, '-append', 'delimiter', '\t');
            clear spikeMat;
        end
        if plot_histogram
            fprintf(histsFid, '\n');
            fclose(histsFid);
            dlmwrite(histsFile, histMat, '-append', 'delimiter', '\t');
            clear histMat;
        end
        if plot_smoothed_histogram
            fprintf(smthhistsFid, '\n');
            fclose(smthhistsFid);
            dlmwrite(smthhistsFile, smthhistMat, '-append', 'delimiter', '\t');
            clear smthhistMat
        end
        if plot_spike_rates
            fprintf(spikeRateFid, '\n');
            fclose(spikeRateFid);
            dlmwrite(spikeRateFile, spikeRateMat, '-append', 'delimiter', '\t');
            clear spikeRateMat;
        end
        if plot_stimulus
           fprintf(stimFid, '\n');
           fclose(stimFid);
           dlmwrite(stimFile, stimMat,  '-append', 'delimiter', '\t');
           clear stimMat;
        end
        if plot_dB_SPL
            fprintf(dBFid, '\n');
            fclose(dBFid);
            dlmwrite(dBFile, dBMat, '-append', 'delimiter', '\t');
            clear dBMat;
        end
    end
    
    function [stimIdx stimPlot] = processStimulus(test_num, trace_num)
        traceData = experiment_data.test(test_num).trace(trace_num);
        if plot_histogram || plot_smoothed_histogram || plot_spike_rates 
            max_hist = max(max(histogram),1);
            offset = max_hist + (0.1 * max_hist);
            for stim_num = 1:size(traceData.stimulus,2)
                stimulus_begin = traceData.stimulus(stim_num).delay;
                if ~isempty(individual_stimulus_signals{test_num,trace_num,stim_num})
                    stim4plot = individual_stimulus_signals{test_num,trace_num,stim_num};
                    stimulus_end = stimulus_begin + ...
                                   length(stim4plot)/stimulus_sampling_frequencies{test_num,trace_num}*1000;
                    stimulus_end = min(stimulus_end,traceData.record_duration);
                    plot_idx = linspace(stimulus_begin,stimulus_end,size(stim4plot,2));
                    %Normalize stimulus to 1
                    stim4plot = stim4plot / max(stim4plot);
                    %Shrink range of stimulus for scale of histogram
                    stim4plot = stim4plot * 0.04 * max_hist;
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    stimIdx = plot_idx;
                    stimPlot = stim4plot;
                else
                    stimulus_end = stimulus_begin + traceData.stimulus(stim_num).duration -1;
                    stim4plot = zeros(1,traceData.stimulus(stim_num).duration);
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    stimIdx = stimulus_begin:stimulus_end;
                    stimPlot = stim4plot;                    
                end
                offset = offset + (0.1 * max_hist);
            end
        else
            offset = 0;
            for stim_num = 1:size(traceData.stimulus,2)
                stimulus_begin = traceData.stimulus(stim_num).delay;
                stimulus_end = stimulus_begin + traceData.stimulus(stim_num).duration - 1;
                if ~isempty(individual_stimulus_signals{test_num,trace_num,stim_num})
                    stim4plot = individual_stimulus_signals{test_num,trace_num,stim_num};
                    plot_idx = linspace(stimulus_begin,stimulus_end,size(stim4plot,2));
                    %Normalize stimulus to 0.9
                    stim4plot = (stim4plot / max(stim4plot))*0.9;
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    stimIdx = plot_idx;
                    stimPlot = stim4plot;
                else
                    stimulus_end = stimulus_begin + traceData.stimulus(stim_num).duration -1;
                    stim4plot = zeros(1,traceData.stimulus(stim_num).duration);
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    stimIdx = stimulus_begin:stimulus_end;
                    stimPlot = stim4plot;      
                end
                offset = offset + 2;
            end
        end
    end
end