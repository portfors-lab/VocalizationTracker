function [individual_stimulus_signals ...
          summed_stimulus_signals ...
          stimulus_sampling_frequencies] = GenerateStimulus(experiment_data, ...
                                                            prefs, ...
                                                            test_nums, ...
                                                            trace_nums)
%
%function [individual_stimulus_signals
%          summed_stimulus_signals
%          stimulus_sampling_frequencies] = GenerateStimulus(experiment_data,
%                                                            prefs,
%                                                            test_nums,
%                                                            trace_nums)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   test_nums           The numbers of the tests to generate the stimulus
%                       for.
%   trace_nums          The numbers of the traces to generate the stimulus
%                       for.
%                       Default: all traces.
%
%   OUTPUT ARGUMENTS
%   individual_stimulus_signals     The generated individual stimulus
%                                   signal components in a cell 
%                                   array indexed by {test_num,trace_num}
%   summed_stimulus_signals         The summed stimulus signals in a cell 
%                                   array indexed by {test_num,trace_num}
%   stimulus_sampling_frequencies   The stimulus sampling frequencies in a 
%                                   cell array indexed by {test_num,trace_num}
warn = 'orange';

if exist('trace_nums','var')
    if isempty(trace_nums)
        trace_nums = [];
    end
else
    trace_nums = [];
end

remove_vocalization_mean = false;
add_rise_fall_to_vocalizations = false;
vocalization_rise_fall_time = 5; %ms

individual_stimulus_signals = cell(max(test_nums),max(trace_nums));
summed_stimulus_signals = cell(max(test_nums),max(trace_nums));
stimulus_sampling_frequencies = cell(max(test_nums),max(trace_nums));

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
        if strcmp(test.testtype,'vocalization')
            stimulus = trace.stimulus;
            if ~isempty(stimulus)
                stim = stimulus(1);
                if ~exist([prefs.audio_directory stim.vocal_call_file], 'file')
                    WriteStatus(['Stimulus file ' prefs.audio_directory stim.vocal_call_file ' not found'...
                        '  -Make sure your stimulus location folder is set correctly'], 'red');
                    return
%                     error(['Stimulus file ' prefs.audio_directory stim.vocal_call_file ' not found']);
                end
                [vocalization sampling_frequency] = ParseAudioData([prefs.audio_directory stim.vocal_call_file]);
                if remove_vocalization_mean
                    vocalization = vocalization - mean(vocalization);
                end
                if stim.reverse_vocal_call
                    vocalization = fliplr(vocalization);
                end
                if ~(sampling_frequency == trace.samplerate_da)
                    WriteStatus('Samplerate of vocalization does not match DA samplerate', warn);
                end

                if add_rise_fall_to_vocalizations
                    %Get length of onset period in samples
                    stim.rise_fall = 5;
                    rise_len = round(sampling_frequency * stim.rise_fall / 1000);
                    rise_envelope = (1-cos(linspace(0,pi,rise_len)))/2;
                    vocalization(1:rise_len) = vocalization(1:rise_len) .* rise_envelope;
                    fall_envelope = fliplr(rise_envelope);
                    vocalization(end-rise_len+1:end) = vocalization(end-rise_len+1:end) .* fall_envelope;
                end
                
                stimulus_length = trace.record_duration / 1000; %Convert to seconds
                stimulus_signal = zeros(1,ceil(stimulus_length * sampling_frequency));
                vocal_delay = stim.delay;
                %Calculate length of the vocalization, in milliseconds
%                 vocal_length = length(vocalization) / sampling_frequency * 1000;
%                 experiment_data.test(test_num).trace(trace_num).stimulus(1).duration = vocal_length;
                stim_offset = round(sampling_frequency * (vocal_delay /1000)) + 1;
                room_at_end = length(stimulus_signal) - (length(vocalization) + stim_offset);
                if room_at_end < 0
                    vocalization = vocalization(1:end+room_at_end);
                end
                stim_samples = length(vocalization);
                %Attenuate the signal
                vocalization = dbAttenuate(vocalization, stim.attenuation);
                %Store the vocalization alone in the stimulus structure
                individual_stimulus_signals{test_num,trace_num,1} = vocalization;
                %Add this signal to the combined stimulus signal
                stimulus_signal(stim_offset:(stim_offset+stim_samples-1)) = ...
                                stimulus_signal(stim_offset:(stim_offset+stim_samples-1)) + vocalization;
            else
                sampling_frequency = trace.samplerate_da;
                stimulus_length = trace.record_duration / 1000; %Convert to seconds
                %Make a base signal the length of the recording duration that
                %the tones will be added to.
                stimulus_signal = zeros(1,round(stimulus_length * sampling_frequency));
            end
        elseif strcmp(test.testtype,'tone') || strcmp(test.testtype,'twotone') || strcmp(test.testtype,'control')
            sampling_frequency = trace.samplerate_da;
            stimulus_length = trace.record_duration / 1000; %Convert to seconds
            %Make a base signal the length of the recording duration that
            %the tones will be added to.
            stimulus_signal = zeros(1,round(stimulus_length * sampling_frequency));
            stimulus = trace.stimulus;  
            if ~isempty(stimulus)
                for stim_num = 1:size(stimulus,2)
                    if strcmp(stimulus(stim_num).soundtype_name,'tone') ||...
                       strcmp(stimulus(stim_num).soundtype_name,'twotone')
                        stim = stimulus(stim_num);
                        stim_time = 0:(1/sampling_frequency):(stim.duration / 1000);
                        stim_signal = sin(2 * pi * stim.frequency * stim_time);
                        if stim.rise_fall > 0
                            %Get length of onset period in samples
                            rise_len = round(sampling_frequency * stim.rise_fall / 1000);
                            rise_envelope = (1-cos(linspace(0,pi,rise_len)))/2;
                            stim_signal(1:rise_len) = stim_signal(1:rise_len) .* rise_envelope;
                            fall_envelope = fliplr(rise_envelope);
                            stim_signal(end-rise_len+1:end) = stim_signal(end-rise_len+1:end) .* fall_envelope;
                        end
                        %Attenuate the signal
                        stim_signal = dbAttenuate(stim_signal, stim.attenuation);
                        stim_offset = sampling_frequency * (stim.delay /1000) + 1;
                        room_at_end = length(stimulus_signal) - (length(stim_signal) + stim_offset);
                        if room_at_end < 0
                            stim_signal = stim_signal(1:end+room_at_end);
                        end
                        stim_samples = length(stim_signal);
                        %Add this signal to the combined stimulus signal
                        stimulus_signal(round(stim_offset):round(stim_offset+stim_samples-1)) = ...
                            stimulus_signal(round(stim_offset):round(stim_offset+stim_samples-1)) + stim_signal;
                        %Save the individual signal for plotting purposes.
                        individual_stimulus_signals{test_num,trace_num,stim_num} = stim_signal;
                    else
                        error(['Unsupported stimulus type: ' stimulus(stim_num).soundtype_name]);
                    end
                end
            end
        else
            %error(['Unsupported test type: ' test.testtype]);
            WriteStatus(['Test ' num2str(test_num) ': Unsupported stimulus test type: ' test.testtype], warn);
            continue
        end
        summed_stimulus_signals{test_num,trace_num} = stimulus_signal;
        stimulus_sampling_frequencies{test_num,trace_num} = sampling_frequency;
    end
end