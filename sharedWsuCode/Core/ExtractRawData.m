function raw_data = ExtractRawData(experiment_data, ...
                                   prefs, ...
                                   test_nums, ...
                                   trace_nums)
%
%function raw_data = ExtractRawData(experiment_data,
%                                   prefs,
%                                   test_nums,
%                                   trace_nums)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   test_nums           The number of the tests to extract data for.
%   trace_nums          The number of the traces to extract data for.
%                       Default: all traces
%
%   OUTPUT ARGUMENTS
%   raw_data            The raw data returned in a cell array indexed 
%                       by {test_num,trace_num}

if exist('test_nums','var')
    if isempty(trace_nums)
        test_nums = 1:length(experiment_data.test);
    end
else
    test_nums = 1:length(experiment_data.test);
end
if exist('trace_nums','var')
    if isempty(trace_nums)
        trace_nums = [];
    end
else
    trace_nums = [];
end

% Author specified flags
normalize_raw_data = false;
remove_mean_raw_data = true;

%Open the raw data file
fid = fopen(prefs.raw_data_filepath,'r','l'); %Open the file, change 'l' -> 'b' for bigEndian byteorder

raw_data = cell(max(test_nums),max(trace_nums));

for test_num = test_nums
    test = experiment_data.test(test_num);
    offset_in_raw_file = test.offset_in_raw_file;
    %Rewind the file
    fseek(fid, 0, 'bof');
    %Seek to the offset of the test data we want
    fseek(fid, offset_in_raw_file, 'bof'); %Seek to the offset

    traces = test.trace;
    if isempty(trace_nums) || max(trace_nums) > length(traces)
        traces_2_process = 1:size(traces,2);
    else
        traces_2_process = trace_nums;
    end
    for trace_num = traces_2_process
        trace = traces(trace_num);
        sample_rate = trace.samplerate_ad;
        num_sweeps = trace.num_samples;
        record_duration_per_run = trace.record_duration;
        samples_per_trace = (record_duration_per_run / 1000) * sample_rate;
        %Rewind the file
        fseek(fid, 0, 'bof');
        %Seek to the offset of the test data we want
        fseek(fid, offset_in_raw_file+((trace_num-1)*(samples_per_trace*num_sweeps)*2), 'bof');
        %Read the raw data segment from the monolithic binary file
        trace_data = fread(fid, [samples_per_trace, num_sweeps], 'int16')';
        if remove_mean_raw_data
            trace_data = trace_data - repmat(mean(trace_data')',1,samples_per_trace);
        end
        if normalize_raw_data
            trace_data = trace_data./max(max(abs(trace_data)));
        else
            trace_data = trace_data./2^15;
        end
        raw_data{test_num,trace_num} = trace_data;
        experiment_data.test(test_num).trace(trace_num).raw_data_samples = length(trace_data);
    end
end

fclose(fid);








