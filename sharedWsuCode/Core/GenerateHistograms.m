function [out_histograms ...
          out_sweep_histograms ...
          out_histogram_bin_widths ...
          out_histogram_bin_centers] = GenerateHistograms(experiment_data, ...
                                                          prefs, ...
                                                          test_nums, ...
                                                          trace_nums, ...
                                                          spike_times)
%
%function [out_histograms
%          out_sweep_histograms
%          out_histogram_bin_widths
%          out_histogram_bin_centers] = GenerateHistograms(experiment_data,
%                                                          prefs,
%                                                          test_nums,
%                                                          trace_nums,
%                                                          bin_width,
%                                                          spike_times)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   test_nums           The numbers of the tests to calculate spike times
%                       for.
%   trace_nums          The numbers of the traces to calculate spike times
%                       for.
%                       Default: all traces
%   spike_times         Pre-calculated spike rates, stored in a cell array 
%                       indexed by {test_num,trace_num}. If empty, the
%                       spike rates are extracted on the fly.
%                       Default = []
%
%   OUTPUT ARGUMENTS
%   experiment_data     bat2matlab structure with requested histograms
%                       included

bin_width = prefs.histogram_bin_width;

if ~exist('trace_nums','var')
    trace_nums = [];
end

if exist('spike_times','var')
    if isempty(spike_times)
        [spike_times spike_idx] = CalculateSpikeTimes(experiment_data,prefs,test_nums,trace_nums);
    end
else
    [spike_times spike_idx] = CalculateSpikeTimes(experiment_data,prefs,test_nums,trace_nums);
end

out_histograms = cell(max(test_nums),max(trace_nums));
out_sweep_histograms = cell(max(test_nums),max(trace_nums));
out_histogram_bin_widths = cell(max(test_nums),max(trace_nums));
out_histogram_bin_centers = cell(max(test_nums),max(trace_nums));

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
        num_bins = ceil(trace.record_duration/bin_width); %Bin size given in milliseconds
        num_sweeps = size(spike_times{test_num,trace_num},1);
        bin_centers = linspace(bin_width/2,trace.record_duration-bin_width/2,num_bins);
        collated_spike_times = [];
        sweep_histograms = zeros(num_sweeps,num_bins);
        trace_spike_times = spike_times{test_num,trace_num};
        for sweep_num = 1:num_sweeps
            sweep_histograms(sweep_num,:) = hist(trace_spike_times{sweep_num},bin_centers);
            collated_spike_times = [collated_spike_times ; trace_spike_times{sweep_num}];
        end
        histogram = hist(collated_spike_times,bin_centers);
        
        out_histograms{test_num,trace_num} = histogram;
        out_sweep_histograms{test_num,trace_num} = sweep_histograms;
        out_histogram_bin_widths{test_num,trace_num} = bin_width;
        out_histogram_bin_centers{test_num,trace_num} = bin_centers;
    end
end

