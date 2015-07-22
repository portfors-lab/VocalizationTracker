function fundamental = SpectrogramFrequencyTracker(y,sampleRate,window_length,frequency_mean)

original_length = length(y);

% y = decimate(y,4);
% sampleRate = round(sampleRate/4);

plot_freq_samples = 2^11;
% plot_time_samples = 2^11;
% plot_time_samples = 2^13;
plot_time_samples = length(y);
% window_length = 2^9;
wl = window_length / sampleRate;
[Hf Ha Ht] = TrackHarmonicsNPS(y, ...
                               sampleRate, ...
                               'nTimes',plot_time_samples, ...
                               'nFrequencies',plot_freq_samples, ...
                               'windowLength',wl, ...
                               'spectrogramType',1, ...
                               ...'nHarmonics',9, ...
                               'nHarmonics',1, ...
                               'plotType', 0, ...
                               'frequencyRange',[5000 20000]);

transient_time = wl/2;
harmonic_sample_rate = plot_time_samples/Ht(end);
harmonic_transient_samples = ceil(harmonic_sample_rate * transient_time);

time_idx = linspace(Ht(1),Ht(end),original_length);
% bad_samples = [];
% next_compare = harmonic_transient_samples - 1;
%Track dominant 1st harmonic, which we know is the fundamental
harmonic = Hf(:);
% harmonic = harmonic./2; %Use these lines when tracking > 1 harmonic
harmonic_time = Ht;


% for i = harmonic_transient_samples:length(harmonic)-harmonic_transient_samples
%     tracked_harmonic = round(harmonic(i)/frequency_mean);
%     if tracked_harmonic == 0
%         harmonic(i) = frequency_mean;
%     else
%         harmonic(i) = harmonic(i)/tracked_harmonic;
%     end
% end
% 
% bad_indices = [];
% for i = harmonic_transient_samples:length(harmonic)-harmonic_transient_samples
%     if harmonic(i) > 1.25*harmonic(next_compare) || harmonic(i) < 0.75*harmonic(next_compare)
%         harmonic(i) = nan;
%         harmonic_time(i) = nan;
%         bad_indices = [bad_indices i];
%     else
%         next_compare = i;
%     end
% end
% good_indices = find(~isnan(harmonic));
% harmonic = harmonic(good_indices);
% harmonic_time = harmonic_time(good_indices);

fundamental = interp1(harmonic_time,harmonic,time_idx,'spline'); 


% fundamental = resample(harmonic,4,1);
