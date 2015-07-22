function [spectrogram_intervals ...
          sampled_frequencies ...
          two_tone_excitory_info] = CalculateSpectrogramIntervals(experiment_data, ...
                                                                  prefs, ...
                                                                  test_nums, ...
                                                                  trace_nums, ...
                                                                  max_interval_width, ...
                                                                  extract_2tone_info)
%
%[spectrogram_intervals
% sampled_frequencies
% two_tone_excitory_info] = CalculateSpectrogramIntervals(experiment_data,
%                                                         prefs,
%                                                         test_nums,
%                                                         trace_nums,
%                                                         max_interval_width)
%
%   INPUT ARGUMENTS
%   experiment_data         Bat2Matlab data structure
%   prefs                   Bat2Matlab preferences
%   test_nums               The numbers of the tests to calculate spectrogram
%                           intervals for.
%   trace_nums              The numbers of the traces to calculate spectrogram
%                           intervals for.
%                           Default: all traces
%   max_interval_width      The maximum interval width that can be returned.
%   extract_2tone_info      Indicates whether to extract the information
%                           about the excitiation stimulus used in two tone
%                           tests.
%                           Default: false
%
%   OUTPUT ARGUMENTS
%   spectrogram_intervals   The spectrogram_intervals that will be used
%                           integrating the power spectral density of audio
%                           stimulus.
%   sampled_frequencies     The unique and sorted frequencies used in the tone
%                           or two tone tests. For vocalization tests, the
%                           fixed intervals are defined in prefs.
%   two_tone_excitory_info  If test_nums includes the numbers of one or
%                           more two tone tests, this return argument
%                           includes the frequency and attenuation of the
%                           excitory signal used in the test.

extract_two_tone_excitory_info = true;

if ~exist('max_interval_width','var')
    max_interval_width = [];
end

sampled_frequencies = [];
two_tone_excitory_info = [];
sampled_frequencies_two_tone_excitory = [];
sampled_attenuations_two_tone_excitory = [];
for test_num = test_nums
    test = experiment_data.test(test_num);
    traces = test.trace;
    num_traces = length(traces);
    if isempty(trace_nums)
        traces_2_process = 1:size(traces,2);
    else
        traces_2_process = trace_nums;
    end
    
    trace_frequencies = [];
    if strcmp(test.testtype,'tone')
        testtype = 1;
    elseif strcmp(test.testtype,'twotone')
        if extract_two_tone_excitory_info
            testtype = 2;
            if length(traces) > 1
                %Find frequency and attenuation of excitory signal in twotone tests
                last_trace = traces(num_traces);
                second_trace = traces(2);
                third_trace = traces(3);
                if length(second_trace.stimulus) ~= 1
                    error(['Second trace of two tone test ' int2str(test_num) ' does not have single tone.']);
                    continue;
                end
                if length(third_trace.stimulus) ~= 1
                    error(['Third trace of two tone test ' int2str(test_num) ' does not have single tone.']);
                    continue;
                end
                if length(last_trace.stimulus) ~= 2
                    error(['Last trace of two tone test ' int2str(test_num) ' does not have two tones.']);
                    continue;
                end

                %Locate the excitory signal test
                if last_trace.stimulus(1).frequency == second_trace.stimulus(1).frequency && last_trace.stimulus(1).attenuation == second_trace.stimulus(1).attenuation
                    excitory_frequency = second_trace.stimulus(1).frequency;
                    excitory_attenuation = second_trace.stimulus(1).attenuation;
                    two_tone_changing_idx = 2;
                elseif last_trace.stimulus(2).frequency == third_trace.stimulus(1).frequency && last_trace.stimulus(2).attenuation == third_trace.stimulus(1).attenuation
                    excitory_frequency = third_trace.stimulus(1).frequency;
                    excitory_attenuation = third_trace.stimulus(1).attenuation;
                    two_tone_changing_idx = 1;
                else
                    display(['No excitory signal found for two tone test ' int2str(test_num)]);
                    continue;
                end
                trace_frequencies_two_tone_excitory = excitory_frequency;
                trace_attenuations_two_tone_excitory = excitory_attenuation;
            end
        else
            testtype = 1;
        end
    elseif strcmp(test.testtype,'control')
        testtype = 1;
    elseif strcmp(test.testtype,'vocalization')
        sampled_frequencies = [10000 : 2000 : 100000];
        testtype = 1;
        continue;
    else
         error(['Test ' int2str(test_num) ' has an unsupported test type for frequency interval calculation.']);
    end
    trace_frequencies = [];
    for trace_num = traces_2_process
        trace = traces(trace_num);
        if testtype == 2 && trace_num < 4
            %First trace is spontaneous rate
            %Second trace is excitation tone
            %Third trace is first paired tone, alone.
            continue;
        end
        if ~isempty(trace.stimulus) && ~strcmp(trace.is_control,'True')
            if length(trace.stimulus) == 0
                error('Debug if we see this');
            elseif length(trace.stimulus) == 1
                stim_idx = 1;
                trace_frequencies = [trace_frequencies trace.stimulus(stim_idx).frequency];
            else %Two tone test
                stim_idx = 1;
                trace_frequencies = [trace_frequencies trace.stimulus(stim_idx).frequency];
                stim_idx = 2;
                trace_frequencies = [trace_frequencies trace.stimulus(stim_idx).frequency];
            end
        end
    end
    sampled_frequencies = [sampled_frequencies trace_frequencies];
    if testtype == 2 && extract_two_tone_excitory_info
        sampled_frequencies_two_tone_excitory = [sampled_frequencies_two_tone_excitory trace_frequencies_two_tone_excitory];
        sampled_attenuations_two_tone_excitory = [sampled_attenuations_two_tone_excitory trace_attenuations_two_tone_excitory];
    end
end

%Sort unique values of all sampled frequencies
freqs = sort(unique(sampled_frequencies));
if testtype == 2 && extract_two_tone_excitory_info
    %Make sure there weren't multiple two-tone tests mixed in.
    two_tone_excitory_info = unique([sampled_frequencies_two_tone_excitory ; ...
                                     sampled_attenuations_two_tone_excitory]','rows');
    if size(unique(two_tone_excitory_info(:,1)),1) ~= 1
        display('Multiple excitation frequencies present in two-tone tests');
    end
end

spectrogram_intervals = [];
if ~isempty(freqs)
    spectrogram_intervals = zeros(length(freqs),2);
    %Calculate lower bound of the first interval
    if isempty(max_interval_width)
        spectrogram_intervals(1,1) = prefs.spectrogram_range(1);
    else
        spectrogram_intervals(1,1) = max(freqs(1)-max_interval_width/2,prefs.spectrogram_range(1));
    end
    if length(freqs) > 1
        for freq_num = 2:length(freqs)
            if isempty(max_interval_width)
                spectrogram_intervals(freq_num-1,2) = (freqs(freq_num)+freqs(freq_num-1))/2;
                spectrogram_intervals(freq_num,1) = (freqs(freq_num)+freqs(freq_num-1))/2;
            else
                spectrogram_intervals(freq_num-1,2) = min(freqs(freq_num-1)+max_interval_width/2,...
                                                       (freqs(freq_num)+freqs(freq_num-1))/2);  
                spectrogram_intervals(freq_num,1) = max(freqs(freq_num)-max_interval_width/2,...
                                                       (freqs(freq_num)+freqs(freq_num-1))/2);  
            end
        end
    end
    %Calculate upper bound of the last interval
    if isempty(max_interval_width)
        spectrogram_intervals(length(freqs),2) = prefs.spectrogram_range(2);
    else
        spectrogram_intervals(length(freqs),2) = min(freqs(end)+max_interval_width/2,prefs.spectrogram_range(2));
    end
end

sampled_frequencies = freqs;





