function integrated_spectral_densities = IntegrateSpectralDensity(spectrogram, ...
                                                                  PSD_freq_range, ...
                                                                  frequency_intervals, ...
                                                                  integration_mode)
%
%function integrated_spectral_densities = IntegrateSpectralDensity(spectrogram,
%                                                                  PSD_freq_range,
%                                                                  frequency_intervals)
%   INPUT ARGUMENTS
%   spectrogram             The spectrogram to integrate over
%   PSD_freq_range          The frequency range of the spectrogram
%   frequency_intervals     The intervals to integrate over. May be in two
%                           forms: [1 2 3] and [1 2 ; 2 3] are equivalent
%                           and indicate two intervals [1 2] and [2 3]
%                           Default = PSD_freq_range (integrate whole
%                           frequency range)
%   integration_mode        Mode = 0: Rectangular integration
%                           Mode = 1: Gaussian integration
%                           Default = 0
%
%   OUTPUT ARGUMENTS
%   integrated_spectral_densities   Each (F X 1) PSD (constant time slice of the
%                                   (F X T) spectrogram) is integrated over each
%                                   frequency interval.
%                                   

if exist('integration_mode','var')
    if isempty(integration_mode)
        integration_mode = 0;
    end
else
    integration_mode = 0;
end

% integration_mode = 0; %Rectangular
% integration_mode = 1; %Gaussian

if exist('frequency_intervals','var')
    if isempty(frequency_intervals)
        frequency_intervals = PSD_freq_range;
    end
else
    frequency_intervals = PSD_freq_range;
end

[num_freq_samples, num_time_samples] = size(spectrogram);
all_freqs = linspace(0,PSD_freq_range(2),num_freq_samples);

%Convert frequency intervals from [1 2 3] format to [1 2 ; 2 3], if
%necessary.
if size(frequency_intervals,2) > 2
    frequency_intervals = frequency_intervals(:);
    frequency_intervals = [frequency_intervals(1:(end-1)) frequency_intervals(2:end)];
end

num_freq_intervals = size(frequency_intervals,1);

%Calculate frequency sample ranges corresponding to the specified intervals
integrated_spectral_densities = zeros(num_freq_intervals,num_time_samples);
freq_sample_start = zeros(num_freq_intervals,1);
freq_sample_end = zeros(num_freq_intervals,1);
for interval_num = 1:num_freq_intervals
    if integration_mode == 0
        %Find the sample in the frequency domain that matches the lower bound
        %of the frequency interval we are integrating over.
        if frequency_intervals(interval_num,1) == 0
            freq_sample_start(interval_num) = 1;
        else
            freq_sample_start(interval_num) = find(frequency_intervals(interval_num,1) >= all_freqs,1,'last');
        end
        %Find the sample in the frequency domain that matches the upper bound
        %of the frequency interval we are integrating over.
        if frequency_intervals(interval_num,2) == all_freqs(end)
            freq_sample_end(interval_num) = length(all_freqs);
        else
            freq_sample_end(interval_num) = find(frequency_intervals(interval_num,2) >= all_freqs,1,'last')-1;
        end
        if (freq_sample_end(interval_num) < freq_sample_start(interval_num))
            integrated_spectral_densities(interval_num) = 0;
        else
            freq_samples = freq_sample_start(interval_num):freq_sample_end(interval_num);
            %Integrate over the power densities where
            %(sampling_frequency/2)./(num_freq_samples-1)
            %is equal to dw : the integration interval
            integrated_spectral_densities(interval_num,:) = sum(spectrogram(freq_samples,:)).*(PSD_freq_range(2))./(num_freq_samples-1);
        end
    else
        mean = (frequency_intervals(interval_num,1) + frequency_intervals(interval_num,2))/2;
        stdev = (frequency_intervals(interval_num,2) - frequency_intervals(interval_num,1))/2;
        distribution = gauss(mean,stdev^2,all_freqs');
        distribution = distribution./max(distribution);
        scaled_spectrogram = spectrogram.*repmat(distribution,1,num_time_samples);
        %Multiply by integration interval
        scaled_spectrogram = scaled_spectrogram.*(PSD_freq_range(2))./(num_freq_samples-1);
        integrated_spectral_densities(interval_num,:) = sum(scaled_spectrogram);
    end
end
     
