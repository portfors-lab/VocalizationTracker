function experiment = ParsePST(filepathname)

%Open the file for I/O
fid = fopen(filepathname);
%Scan each line into a cell array
% file = textscan(fid, '%s', 'delimiter', '\n');
file = textscan(fid, '%s', 9,'delimiter', '\n');

%Get the first column of the cell array (the lines)
lines = file{1};
%Close the file for I/O


%Static strings
test_aborted = 'Test aborted.';
end_id = 'End of ID information';
end_test_parameters = 'End of test parameters';
end_spike_data = 'End of spike data';
end_auto_test = 'End of auto test';

test_types = {'tone',...
              'fmsweep',....
              'synthesized_batsound',....
              'amsound',....
              'broad_band_noise',....
              'narrow_band_noise',....
              'click',....
              'vocalization',....
              'high_pass_noise',....
              'low_pass_noise',....
              'sine_wave_modulation',....
              'square_wave_modulation'};

%Indicates the current line being read.
line_num = 1;
%Total lines in the PST file
num_lines = length(lines);
%Offset into the raw data file
raw_pos = 0;

%Collect experiment-wide data from ID section
experiment.pst_filename = lines{1};
experiment.date = lines{2};
experiment.title = lines{3};
experiment.who = lines{4};
experiment.computername = lines{5};
experiment.program_date = lines{6};

%Scan past the ID section
while ~strcmp(lines{line_num},end_id)
    line_num = line_num + 1;
end
line_num = line_num + 1;

%Iterate through each test
line = fgets(fid);
while line ~=-1
    
    %Get the test type
%     full_testtype = lines{line_num};    
    full_testtype=line;
%     line_num = line_num + 1;
    %Get the test number.
    line = fgets(fid);
    test_line = textscan(line,'%n %*s %*s %*n %s',1);
    test_num = test_line{1};
    test_time=test_line{2};
    experiment.test(test_num).testnum = test_num;
    experiment.test(test_num).time=test_time{1};
    %This is the Batlab assigned test type, not the test type that Bat2Matlab uses
    experiment.test(test_num).full_testtype = full_testtype;
    %Store the beginning location of the test

    %Extract the test paramewters
    line = fgets(fid);
    test_data = textscan(line,'%n %s %n %n %n %n %s %n %n %n %n %n %n',1);
    num_traces = test_data{1};

    
    %Not even sure what these test are. Not necessary?
%     vocal_call_io_test_number = textscan(lines{line_num},'%n');
%     line_num = line_num + 1;
%     vocal_call_io_test_number = vocal_call_io_test_number{1};
%     if vocal_call_io_test_number > 0
%         error('Handle vocal io tests');
%     end
    
    %Scan past the test parameter section
    line = fgets(fid);
    while ~strcmp(strtrim(line),end_test_parameters)
       line = fgets(fid);
    end
    line = fgets(fid);
    %Get the position of the beginning of the test in the raw data file
    experiment.test(test_num).offset_in_raw_file = raw_pos;
    for trace_num = 1:num_traces
        %Get the trace data for this test
        trace_data = textscan(line,'%n');
        trace_data = trace_data{1};
        num_sweeps = trace_data(1);
        samplerate_da = trace_data(2);
        samplerate_ad = trace_data(4);
        duration = trace_data(5);
        points = samplerate_ad/1000*duration;
        experiment.test(test_num).trace(trace_num).record_duration = duration;
        experiment.test(test_num).trace(trace_num).samplerate_da = samplerate_da;
        experiment.test(test_num).trace(trace_num).samplerate_ad = samplerate_ad;
        experiment.test(test_num).trace(trace_num).num_samples = num_sweeps;
        %The length of the trace in the raw data file
        trace_raw_data_length = points*num_sweeps*2;
        %Store the trace raw offset and length
        experiment.test(test_num).trace(trace_num).offset_in_raw_file = raw_pos;
        experiment.test(test_num).trace(trace_num).length_in_raw_file = trace_raw_data_length;
        %Increment the position in the raw data file
        raw_pos = raw_pos + trace_raw_data_length;
        
        %Collect the stimulus parameters for all 4 channels
        stimulus_num = 0;
        test_stim_type = 0; %Default
        
        %skip 3 lines
        textscan(fid, '%s', 3,'delimiter', '\n');
        line=fgets(fid);
        for channel_num = 1:4
%             [stimulus stim_type] = ParsePSTStimulus(lines{line_num - 1 + 5*channel_num},test_num,trace_num);
            [stimulus stim_type] = ParsePSTStimulus(line,test_num,trace_num);
            if ~isempty(stimulus) %&& stim_type ~= 5 
                stimulus_num = stimulus_num + 1;
                experiment.test(test_num).trace(trace_num).stimulus(stimulus_num) = stimulus;
                if test_stim_type == 0
                    test_stim_type = stim_type;
                end
            end
            %skip 5 lines
            textscan(fid, '%s', 4,'delimiter', '\n');
            line=fgets(fid);
        end

        %If there is no stimulus
        if ~isfield(experiment.test(test_num).trace(trace_num), 'stimulus')
            experiment.test(test_num).trace(trace_num).stimulus = [];
        end
        
        %Set the test type based on the stimulus variety
        switch test_stim_type
            case 0, %No Stimulus
                test_type = 'control';
            case 1, %tone
                switch stimulus_num
                    case 1, test_type = 'tone';
                    case 2, test_type = 'twotone';
                    case 3, test_type = 'threetone';
                    case 4, test_type = 'fourtone';
                end
            otherwise, test_type = test_types{test_stim_type};
        end
        experiment.test(test_num).testtype = test_type;
        if test_stim_type > 0
            experiment.test(test_num).trace(trace_num).is_control = 0;
        else
            %If there is no stimulus present in the trace, then label it as a contol.
            experiment.test(test_num).trace(trace_num).is_control = 1;
        end
        
        %Scan past the spike data
        while ~strcmp(strtrim(line),end_spike_data)
           line = fgets(fid);
        end
        line = fgets(fid);
    end
    
    %Assert that the spike data terminates with an end auto test statement
    if ~strcmp(strtrim(line),end_auto_test)
        error('No end of auto test found after spike data');
    end
    
    %Get the test comment
    line = fgets(fid);
    experiment.test(test_num).comment = line;
    
    line = fgets(fid);
end
fclose(fid);


    
    
    
    
    
    
    
    
    
