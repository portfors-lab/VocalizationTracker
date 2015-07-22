function signal_peaks = EnforceRefractoryPeriod(signal_peaks, refractory_period, fs)
    % Function to implement a refractory period on the calculated spike times
    spike_times = (signal_peaks-1)/fs*1000; %in milliseconds
    previous_spike_time = spike_times(1);
    good_spike_indices = zeros(1,length(spike_times));
    good_spike_indices(1) = 1;
    num_good_spikes = 1;
    for i = 2:length(spike_times)
        this_spike_time = spike_times(i);
        if (this_spike_time - previous_spike_time) >= refractory_period
            num_good_spikes = num_good_spikes + 1;
            good_spike_indices(num_good_spikes) = i;
            previous_spike_time = this_spike_time;
        end
    end
    signal_peaks = signal_peaks(good_spike_indices(1:num_good_spikes));
end