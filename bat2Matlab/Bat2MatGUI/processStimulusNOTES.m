function processStimulus

traceData = prefs.test(test).trace(trace);
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
                    stimMat(:,index) = plot_idx;
                    stimMat(:, index+1) = stim4plot;
                else
                    stimulus_end = stimulus_begin + traceData.stimulus(stim_num).duration -1;
                    stim4plot = zeros(1,traceData.stimulus(stim_num).duration);
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    stimMat(:,index) = stimulus_begin:stimulus_end;
                    stimMat(:, index+1) = stim4plot;                    
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
                    stimMat(:,index) = plot_idx;
                    stimMat(:, index+1) = stim4plot;
                else
                    stimulus_end = stimulus_begin + traceData.stimulus(stim_num).duration -1;
                    stim4plot = zeros(1,traceData.stimulus(stim_num).duration);
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    stimMat(:,index) = stimulus_begin:stimulus_end;
                    stimMat(:, index+1) = stim4plot;
                end
                offset = offset + 2;
            end
        end