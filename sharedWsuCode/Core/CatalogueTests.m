function catalogue = CatalogueTests(experiment_data,test_nums,cell_id,spike_threshold)
%
%function catalogue = CatalogueTests(experiment_data,test_nums)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab structure
%   test_nums           Test numbers to catalogue
%
%   OUTPUT ARGUMENTS
%   catalogue           A structure cataloguing pertinent data
%                       from the requested tests. This data is 
%                       redundant with experiment_data, but is
%                       more searchable in this form.

total_vocalization_traces = 0;
total_tone_traces = 0;
total_twotone_traces = 0;
total_traces = 0;
for test_num = test_nums
    test = experiment_data.test(test_num);
    total_traces = total_traces + length(size(test.trace,2));
    switch test.testtype
        case 'vocalization'
            total_vocalization_traces = total_vocalization_traces + length(size(test.trace,2));
        case 'tone'
            total_tone_traces = total_tone_traces + length(size(test.trace,2));
        case 'twotone'
            total_twotone_traces = total_twotone_traces + length(size(test.trace,2));
    end
end

catalogue = struct('cell_name', nan, ...
                   'test_type', nan, ...
                   'test_num', nan, ...
                   'trace_num', nan, ...
                   'stim1_length', nan, ...
                   'stim1_freq', nan, ...
                   'stim1_attn', nan, ...
                   'stim1_delay', nan, ...
                   'stim2_length', nan, ...
                   'stim2_freq', nan, ...
                   'stim2_attn', nan, ...
                   'stim2_delay', nan, ...
                   'vocal_file_name', nan, ...
                   'record_duration', nan, ...
                   'spike_threshold', nan);
              

trace_idx = 0;

for test_num = test_nums
    test = experiment_data.test(test_num);
    traces = test.trace;
    traces_2_process = 1:size(traces,2);
    for trace_num = traces_2_process
        trace_idx = trace_idx + 1;
        
        trace = traces(trace_num);
        stimulus = trace.stimulus;      
        
        catalogue(trace_idx).cell_name = cell_id;
        catalogue(trace_idx).test_num = test_num;
        catalogue(trace_idx).trace_num = trace_num;
        catalogue(trace_idx).test_type = '';
        catalogue(trace_idx).vocal_file_name = '';
        catalogue(trace_idx).record_duration = trace.record_duration;
        catalogue(trace_idx).spike_threshold = spike_threshold;
        
        if ~isempty(stimulus)
            switch test.testtype
                case 'vocalization'
                    catalogue(trace_idx).test_type = 'vocal';
                    stim = stimulus(1);
                    catalogue(trace_idx).vocal_file_name = lower(stim.vocal_call_file);
                    catalogue(trace_idx).stim1_attn = stim.attenuation;
                    catalogue(trace_idx).stim1_delay = stim.delay;
                case {'tone', 'twotone'}
                    num_stim = size(stimulus,2);
                    if num_stim == 1
                        catalogue(trace_idx).test_type = 'one_tone';
                        stim = stimulus(1);
                        catalogue(trace_idx).stim1_length = stim.duration;
                        catalogue(trace_idx).stim1_freq = stim.frequency;
                        catalogue(trace_idx).stim1_attn = stim.attenuation;
                        catalogue(trace_idx).stim1_delay = stim.delay;
                    elseif num_stim == 2
                        catalogue(trace_idx).test_type = 'two_tone';
                        stim = stimulus(1);
                        catalogue(trace_idx).stim1_length = stim.duration;
                        catalogue(trace_idx).stim1_freq = stim.frequency;
                        catalogue(trace_idx).stim1_attn = stim.attenuation;
                        catalogue(trace_idx).stim1_delay = stim.delay;
                        
                        stim = stimulus(2);
                        catalogue(trace_idx).stim2_length = stim.duration;
                        catalogue(trace_idx).stim2_freq = stim.frequency;
                        catalogue(trace_idx).stim2_attn = stim.attenuation;
                        catalogue(trace_idx).stim2_delay = stim.delay;
                    end
                otherwise
                    error(['Unknown test type: ' test.testtype]);
            end
        end
    end
end

if trace_idx == 0
    catalogue = [];
end


    