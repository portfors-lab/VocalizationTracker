function [out_spike_rates, ...
          out_spike_rate_averages, ...
          out_spike_rate_sampling_frequencies] = GenerateSpikeRates(experiment_data, ...
                                                                    prefs, ...
                                                                    test_nums, ...
                                                                    trace_nums, ...
                                                                    spike_times)
%
%function [out_spike_rates,
%          out_spike_rate_averages,
%          out_spike_rate_sampling_frequencies] = GenerateSpikeRates(experiment_data,
%                                                                    test_nums,
%                                                                    trace_nums,
%                                                                    sampling_frequency,
%                                                                    spike_times)
%
%   INPUT ARGUMENTS
%   experiment_data         Bat2Matlab data structure
%   prefs                   Bat2Matlab preferences
%   test_nums               The numbers of the tests to calculate spike rates
%                           for.
%   trace_nums              The numbers of the traces to calculate spike rates
%                           for.
%                           Default: all traces
%   sampling_frequency      Sampling frequency for rate calculation
%                           Default: Sampling rate of MER recording
%   spike_times             Pre-calculated spike rates, stored in a cell array 
%                           indexed by {test_num,trace_num}. If empty, the
%                           spike rates are extracted on the fly.
%                           Default = []
%
%   OUTPUT ARGUMENTS      
%   out_spike_rates                     The calculated spike times returned in a cell array 
%                                       indexed by {test_num,trace_num}
%   out_spike_rate_averages             The calculated spike rate averagest returned in a cell 
%                                       array indexed by {test_num,trace_num}
%   out_spike_rate_sampling_frequencies The calculated spike rate sampling frequencies 
%                                       returned in a cell array indexed by {test_num,trace_num}

%last modified 5/26/11 Amy Boyle, edited command line output

warnColor = 'orange';

if ~exist('trace_nums','var')
    trace_nums = [];
end

if ~exist('sampling_frequency','var')
    sampling_frequency = [];
end

if exist('spike_times','var')
    if isempty(spike_times)
        spike_times = CalculateSpikeTimes(experiment_data, ...
                                          prefs, ...
                                          test_nums, ...
                                          trace_nums);
    end
else
    spike_times = CalculateSpikeTimes(experiment_data, ...
                                      prefs, ...
                                      test_nums, ...
                                      trace_nums);
end

out_spike_rates = cell(max(test_nums),max(trace_nums));
out_spike_rate_averages = cell(max(test_nums),max(trace_nums));
out_spike_rate_sampling_frequencies = cell(max(test_nums),max(trace_nums));

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

        multSpikesCount=0;
        
        if isnan(prefs.spike_rate_cutoff_frequency)
            sampling_frequency = trace.samplerate_ad;
        else
            sampling_frequency = prefs.spike_rate_sampling_frequency;
        end
        
        bin_width = 1000 / sampling_frequency;
        num_bins = trace.record_duration/bin_width; %Bin size given in milliseconds
        num_sweeps = size(spike_times{test_num,trace_num},1);
        bin_centers = linspace(bin_width/2,trace.record_duration-bin_width/2,num_bins);     
        spike_rates = zeros(num_sweeps, num_bins);
        trace_spike_times = spike_times{test_num,trace_num};
        for sweep_num = 1:num_sweeps
            sweep_histogram = hist(trace_spike_times{sweep_num},bin_centers);
            if max(sweep_histogram) > 1
                %warning('Multiple spikes within time bin for binary spike train calculation');
                multSpikesCount = multSpikesCount+1;
                sweep_histogram(sweep_histogram > 1) = 1;
            end
            spike_rate = RateFilter(sweep_histogram,sampling_frequency,prefs.spike_rate_cutoff_frequency)';
            spike_rates(sweep_num,:) = spike_rate;
        end
        if multSpikesCount>0
            WriteStatus(['Test: ' num2str(test_num) ', Trace: ' num2str(trace_num) ' Multiple spikes within time bin for binary spike train calculation (x' num2str(multSpikesCount)], warnColor);
        end
        out_spike_rates{test_num,trace_num} = spike_rates;
        out_spike_rate_averages{test_num,trace_num} =  mean(spike_rates);
        out_spike_rate_sampling_frequencies{test_num,trace_num} = sampling_frequency;
    end
end

