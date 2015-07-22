function [out_spike_times ...
          out_spike_idxs ...
          out_smoothed_data ...
          out_peak_thresholds] = CalculateSpikeTimes(experiment_data, ...
                                                     prefs, ....
                                                     test_nums, ...
                                                     trace_nums, ...
                                                     raw_data)
%
%[out_spike_times
% out_spike_idxs
% out_smoothed_data
% out_peak_thresholds] = CalculateSpikeTimes(experiment_data,
%                                            prefs,
%                                            test_nums,
%                                            trace_nums,
%                                            peak_threshold,
%                                            filter_cutoff,
%                                            power_exponent,
%                                            raw_data)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   test_nums           The numbers of the tests to calculate spike times
%                       for.
%   trace_nums          The numbers of the traces to calculate spike times
%                       for.
%                       Default: all traces
%   raw_data            Extracted raw data, stored in a cell array 
%                       indexed by {test_num,trace_num}. If empty, the raw
%                       data is extracted on the fly.
%                       Default = []
%
%   OUTPUT ARGUMENTS
%   out_spike_times     The calculated spike times returned in a cell array 
%                       indexed by {test_num,trace_num}
%   out_spike_idxs      The calculated spike time indexes returned in a cell array 
%                       indexed by {test_num,trace_num}
%   out_smoothed_data   The calculated smoothed MER data returned in a cell array 
%                       indexed by {test_num,trace_num}
%   out_peak_thresholds The peak threshold level used returned in a cell array 
%                       indexed by {test_num,trace_num}


filter_cutoff = prefs.spike_time_filter_cutoff;
power_exponent = prefs.spike_time_power_exponent;
peak_threshold = prefs.spike_time_peak_threshold;
refractory_period = prefs.spike_time_refractory_period;

if exist('trace_nums','var')
    if isempty(trace_nums)
        trace_nums = [];
    end
else
    trace_nums = [];
end

out_spike_times = cell(max(test_nums),max(trace_nums));
out_spike_idxs = cell(max(test_nums),max(trace_nums));
out_smoothed_data = cell(max(test_nums),max(trace_nums));
out_peak_thresholds = cell(max(test_nums),max(trace_nums));

for test_num = test_nums
    test = experiment_data.test(test_num);
    traces = test.trace;
    if isempty(trace_nums) || max(trace_nums) > length(traces)
        traces_2_process = 1:size(traces,2);
    else
        traces_2_process = trace_nums;
    end
    for trace_num = traces_2_process
        trace = traces(trace_num);
        fs = trace.samplerate_ad;
        
        if exist('raw_data','var')
            if isempty(raw_data)
                raw_data = ExtractRawData(experiment_data, raw_data_filepath, test_num, trace_num);
                trace_data = raw_data{test_num,trace_num};
                clear raw_data;
            else
                trace_data = raw_data{test_num,trace_num};
            end
        else
            raw_data = ExtractRawData(experiment_data, prefs, test_num, trace_num);
            trace_data = raw_data{test_num,trace_num};
            clear raw_data;
        end
        %[num_sweeps sweep_length] = size(raw_data{test_num, trace_num});
        [num_sweeps sweep_length] = size(trace_data);
        spike_idx = cell(num_sweeps,1);
        spike_times = cell(num_sweeps,1);
        smoothed_data = zeros(num_sweeps, sweep_length);
        for sweep_num = 1:num_sweeps
            signal = trace_data(sweep_num,:);
            %signal = raw_data{test_num, trace_num}(sweep_num,:);
            smoothed_signal = Lowpass(abs(signal).^power_exponent,fs,filter_cutoff,4);  % Smooth the signal.^pw
            smoothed_signal = abs(smoothed_signal).^(1/power_exponent);                 % Scale down the power of signal
            signal_peaks    = DetectMaxima(smoothed_signal,peak_threshold,0);           % Detect maxima
            % Test to see whether the microphone interfered with the electrode recording by finding the ratio
            % of the mean power of the signal before the stimulus delay to the mean power of the signal for a 
            % period of 10ms following the stimulus onset. This is roughly a period before the cell would start
            % responding to the audio stimulus, but where the electrode would be interfered with.
            if ~trace.is_control && length(signal_peaks) > 40
                delay = trace.stimulus.delay;
                start_sim = floor(delay*fs/1000);
                mean_power_before_stim = mean(smoothed_signal(1:start_sim));
                mean_power_after_stim = mean(smoothed_signal(start_sim:start_sim+floor(10*fs/1000)));
                if mean_power_after_stim/mean_power_before_stim > 25
                    warning(['Saturated MER signal detected for test ' int2str(test_num) ' trace ' int2str(trace_num) ' sweep ' int2str(sweep_num)]);
                    signal_peaks = [];
                end
            end
            
            if (refractory_period > 0 && length(signal_peaks) > 0)
                signal_peaks = EnforceRefractoryPeriod(signal_peaks, refractory_period, fs);
            end
            signal_length = length(signal)/fs*1000; % in milliseconds
            max_spikes = signal_length / refractory_period;

            spike_idx{sweep_num,1} = signal_peaks;
            spike_times{sweep_num,1} = (signal_peaks-1)/fs*1000; %in milliseconds
            
            if nargout >2
                %Only cache this data if it is to be returned
                smoothed_data(sweep_num,:) = smoothed_signal;
            end
        end

        out_spike_times{test_num,trace_num} = spike_times;
        out_spike_idxs{test_num,trace_num} = spike_idx;
        if nargout > 2
            out_smoothed_data{test_num,trace_num} = smoothed_data;
            out_peak_thresholds{test_num,trace_num} = peak_threshold;
        end
    end
    
end



