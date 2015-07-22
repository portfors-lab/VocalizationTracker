function [out_stimulus_spectrogram_harmonic_vals ...
          out_stimulus_spectrogram_harmonic_idx] = FindHarmonics(experiment_data, ...
                                                                 prefs, ...
                                                                 test_nums, ...
                                                                 trace_nums, ...
                                                                 stimulus_spectrograms, ...
                                                                 stimulus_spectrogram_freq_idx)
%
%function [out_stimulus_spectrogram_harmonic_vals
%          out_stimulus_spectrogram_harmonic_idx] = FindHarmonics(experiment_data,
%                                                                 prefs,
%                                                                 test_nums,
%                                                                 trace_nums,
%                                                                 stimulus_spectrograms,
%                                                                 stimulus_spectrogram_freq_idx)
%
%   INPUT ARGUMENTS
%   experiment_data                 Bat2Matlab data structure
%   prefs                           Bat2Matlab preferences
%   test_nums                       The numbers of the tests to find the harmonics for.
%   trace_nums                      The numbers of the traces to find the harmonics for.
%                                   Default: all traces.
%   stimulus_spectrograms           The spectrograms to be analyzed, in a
%                                   cell array indexed by {test_num,trace_num}
%                                   If empty, is calculated on the fly
%                                   Default = []
%   stimulus_spectrogram_freq_idx   The spectrograms' frequency idices, in a
%                                   cell array indexed by {test_num,trace_num}
%                                   If empty, is calculated on the fly
%                                   Default = []
%
%   OUTPUT ARGUMENTS
%   out_stimulus_spectrogram_harmonic_vals  The harmonic frequencies returned in a cell 
%                                           array indexed by {test_num,trace_num}
%   out_stimulus_spectrogram_harmonic_idx   The harmonic indices returned in a cell 
%                                           array indexed by {test_num,trace_num}

if exist('trace_nums','var')
    if isempty(trace_nums)
        trace_nums = [];
    end
else
    trace_nums = [];
end

if exist('stimulus_spectrograms','var')
    if isempty(stimulus_spectrograms)
        [stimulus_spectrograms ...
         stimulus_spectrogram_freq_idx] = ExtractRawData(experiment_data, ...
                                                         raw_data_filepath, ...
                                                         test_nums, ...
                                                         trace_nums);
    end
else
    [stimulus_spectrograms ...
     stimulus_spectrogram_freq_idx] = ExtractRawData(experiment_data, ...
                                                     raw_data_filepath, ...
                                                     test_nums, ...
                                                     trace_nums);
end

out_stimulus_spectrogram_harmonic_vals = cell(max(test_nums),max(trace_nums));
out_stimulus_spectrogram_harmonic_idx = cell(max(test_nums),max(trace_nums));

for test_num = test_nums
    test = experiment_data.test(test_num);
    traces = test.trace;
    if isempty(trace_nums)
        traces_2_process = 1:size(traces,2);
    else
        traces_2_process = trace_nums;
    end
    for trace_num = traces_2_process
        trace = traces(trace_num);
        spectrogram = stimulus_spectrograms{test_num,trace_num};
        dB_ref = max(max(spectrogram));
        %Convert to dB scale for calculation, using the maximum value of the
        %spectrogram as the reference.
        warning off
        spectrogram = 10*log10(spectrogram/dB_ref);
        warning on
        dbMax = max(max(spectrogram));
        if dbMax == -Inf
            spectrogram = spectrogram*0;
        else
            %Confine range of spectrogram to dbRange decibels.
            spectrogram = max(spectrogram,-prefs.dbRange);
        end
        
        sample_rate = length(stimulus_spectrogram_freq_idx{test_num,trace_num});
        cutoff = sample_rate/prefs.frequency_cutoff_ratio;
        spectrogram_minimum = min(min(spectrogram));
        %Make spectrogram positive for maxima detection
        spectrogram = spectrogram - spectrogram_minimum;
        %Look for peaks above the noise floor
        min_val = -spectrogram_minimum * prefs.harmonic_peak_detection_noise_floor;
        num_time_samples = size(spectrogram,2);
        num_freq_samples = size(spectrogram,1);
        peak_indices = zeros(prefs.num_harmonics,num_freq_samples);
        peak_vals = zeros(prefs.num_harmonics,num_freq_samples);
        for time_sample_num = 1:num_time_samples
            if time_sample_num == 400
                foo = 1;
            end
            signal = spectrogram(:,time_sample_num);
            if var(signal) == 0
                max_idx = [];
            else
                smoothed_signal = Lowpass(signal,sample_rate,cutoff,4);
                max_idx = DetectMaxima(smoothed_signal,min_val,0);
            end
            num_peaks = length(max_idx);
            max_vals = signal(max_idx);
            
            peaks = zeros(prefs.num_harmonics,2);
            if num_peaks > 0
                peaks2sort = [max_vals(:) max_idx(:)];
                %Sort peaks by amplitude of peaks, in decending order
                sorted_peaks = sortrows(peaks2sort,-1);
                %Prune if there are too many found
                sorted_peaks = sorted_peaks(1:min(num_peaks,prefs.num_harmonics),:);
                %Sort peaks by index
                sorted_peaks = sortrows(sorted_peaks,2);
                peaks(1:min(num_peaks,prefs.num_harmonics),:) = sorted_peaks;
            end
            %Add to signal length matrix to form one signal per hermonic
            peak_vals(:,time_sample_num) = peaks(:,1);
            peak_indices(:,time_sample_num) = peaks(:,2);
        end
        %Save the frequencies and amplitudes
        out_stimulus_spectrogram_harmonic_vals{test_num,trace_num} = peak_vals;
        out_stimulus_spectrogram_harmonic_idx{test_num,trace_num} = peak_indices;
    end
end













