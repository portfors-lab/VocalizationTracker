function [out_stimulus_spectrograms ...
          out_stimulus_spectrogram_time_idxs ...
          out_stimulus_spectrogram_freq_idxs ...
          out_stimulus_dB_SPLs] = GenerateSpectrograms(experiment_data, ...
                                                       prefs, ...
                                                       test_nums, ...
                                                       trace_nums, ...
                                                       stimulus_signals)
%
%function [out_stimulus_spectrograms
%          out_stimulus_spectrogram_time_idxs
%          out_stimulus_spectrogram_freq_idxs
%          out_stimulus_dB_SPLs] = GenerateSpectrograms(experiment_data,
%                                                       prefs,
%                                                       test_nums,
%                                                       trace_nums,
%                                                       stimulus_signals)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   test_nums           The numbers of the tests to generate the
%                       spectrograms for.
%   trace_nums          The numbers of the tests to generate the
%                       spectrograms for.
%                       Default: all traces.
%   stimulus_signals    Pre-calculated stimulus signals, stored in a cell 
%                       array indexed by {test_num,trace_num}. If empty, 
%                       the signals are extracted on the fly.
%                       Default = []
%
%   OUTPUT ARGUMENTS
%   out_stimulus_spectrograms           The generated spectrograms in a cell 
%                                       array indexed by {test_num,trace_num}
%   out_stimulus_spectrogram_time_idxs  The generated spectrogram time indices 
%                                       in a cell array indexed by
%                                       {test_num,trace_num}
%   out_stimulus_spectrogram_freq_idxs  The generated spectrogram frequency 
%                                       indices in a cell array indexed by
%                                       {test_num,trace_num}
%   out_stimulus_dB_SPLs                The stimulus signal on a dB SPL scale 
%                                       in a cell array indexed by
%                                       {test_num,trace_num}

%Modified Amy Boyle March 1 2011, included colormap argument for
%NonparametricSpectrogram

if exist('trace_nums','var')
    if isempty(trace_nums)
        trace_nums = [];
    end
else
    trace_nums = [];
end

if exist('stimulus_signals','var')
    if isempty(trace_nums)
        [individual_stimulus_signals ...
         stimulus_signals] = GenerateStimulus(experiment_data, ...
                                              prefs, ...
                                              test_nums, ...
                                              trace_nums);
    end
else
    [individual_stimulus_signals ...
     stimulus_signals] = GenerateStimulus(experiment_data, ...
                                          prefs, ...
                                          test_nums, ...
                                          trace_nums);
end
clear individual_stimulus_signals;

out_stimulus_spectrograms = cell(max(test_nums),max(trace_nums));
out_stimulus_spectrogram_time_idxs = cell(max(test_nums),max(trace_nums));
out_stimulus_spectrogram_freq_idxs = cell(max(test_nums),max(trace_nums));
out_stimulus_dB_SPLs = cell(max(test_nums),max(trace_nums));

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
        stimulus_signal = stimulus_signals{test_num,trace_num};
        sampling_frequency = trace.samplerate_da;
        stimulus_length = trace.record_duration; %In milliseconds
        spectrogram_time_samples = stimulus_length * prefs.spectrogram_time_samples_per_millisecond;
        if var(stimulus_signal) ~= 0
            [S,times,frequencies] = NonparametricSpectrogram(stimulus_signal,...
                                                             sampling_frequency,...
                                                             'nFrequencies',prefs.spectrogram_freq_samples,...
                                                             'nTimes',spectrogram_time_samples,...
                                                             'windowLength',prefs.spectrogram_window_length,...
                                                             'frequencyRange',prefs.spectrogram_range,...
                                                             'colormap', prefs.colormap_name);
                                                         
        else
            S = zeros(prefs.spectrogram_freq_samples+1,spectrogram_time_samples);
            %Generate indices for plotting axes
            times = linspace(0,stimulus_length,spectrogram_time_samples);
            frequencies = linspace(0,sampling_frequency/2,prefs.spectrogram_freq_samples+1);
        end
        %Convert to power spectral density
        S = abs(S).^2;       
%         load(['speaker_model_' int2str(round(sampling_frequency/1000))],'model');
%         voltage_to_pressure_conversion_coefficients = feval(model,frequencies);
%         S = S.*repmat(voltage_to_pressure_conversion_coefficients(:),1,length(times));
        
        
        rms_pressure_power = IntegrateSpectralDensity(S,prefs.spectrogram_range);
        warning off %Suppress log(0) warnings
        dB_SPL = 10*log10(rms_pressure_power/prefs.dB_SPL_ref);
        warning on
        dB_SPL = max(dB_SPL,prefs.dbMin);
        
        %Save the data for output
        out_stimulus_spectrograms{test_num,trace_num} = S;
        out_stimulus_spectrogram_time_idxs{test_num,trace_num} = times;
        out_stimulus_spectrogram_freq_idxs{test_num,trace_num} = frequencies;
        out_stimulus_dB_SPLs{test_num,trace_num} = dB_SPL;
    end
end













