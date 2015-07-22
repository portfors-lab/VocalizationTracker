close all
clear variables

debug = false;
plot_interpolation = true;

prefs = GeneratePreferences();

%Load the dB SPL callibration file for the speaker
eval(prefs.speaker_calibration_file);
%The different sampling frequencies we need to calibrate for, in Hz
calibration_sampling_frequencies = [250 333 400]*1000;

signal_length = 10;%ms
spectrogram_time_samples = signal_length * prefs.spectrogram_time_samples_per_millisecond;

num_calibrations = length(calibration_sampling_frequencies);
num_freqs = length(sampled_frequencies);
calibration_idx = 0;
calibration_coefficients = zeros(num_freqs,num_calibrations);
for calibration_sampling_frequency = calibration_sampling_frequencies
    calibration_idx = calibration_idx +1;
    freq_idx = 0;
    %Cycle through the frequencies sampled in the speaker calibration
    for sampled_frequency = sampled_frequencies
        freq_idx = freq_idx + 1;
        signal_amplitude = 1;
        time_idx = 0:1/calibration_sampling_frequency:signal_length/1000;
        %Generate a signal
        signal = signal_amplitude*sin(2*pi*sampled_frequency*time_idx);
        [S,t,f] = NonparametricSpectrogram(signal,...
                                           calibration_sampling_frequency,...
                                           'nFrequencies',prefs.spectrogram_freq_samples,...
                                           'nTimes',spectrogram_time_samples,...
                                           'windowLength',prefs.spectrogram_window_length,...
                                           'plotType',0,...
                                           'frequencyRange',prefs.spectrogram_range);

        time_slice = abs(S(:,round(spectrogram_time_samples/2)));
        %Convert to power spectral density
        time_slice = time_slice.^2;
        %Calculate rms pressure power over the time slice
        time_slice_rms = IntegrateSpectralDensity(time_slice,prefs.spectrogram_range);
        %Calculate the conversion coefficient required to achieve desired
        %loudness in dB SPL
        volt_to_pressure_conversion = exp(log(10)*(max_dB_SPL(freq_idx)/10-log10(time_slice_rms)+log10(prefs.dB_SPL_ref)));
        calibration_coefficients(freq_idx,calibration_idx) = volt_to_pressure_conversion;
        if debug
            sampled_frequency
            time_slice = time_slice * volt_to_pressure_conversion;
            time_slice_in_dBs = 10*log10(time_slice/prefs.dB_SPL_ref);
            max_of_time_slice = max(time_slice_in_dBs)
            time_slice_rms = IntegrateSpectralDensity(time_slice,prefs.spectrogram_range);
            target_dB_SPL = max_dB_SPL(freq_idx)
            calculated_dB_SPL = 10*log10(time_slice_rms/prefs.dB_SPL_ref)
        end
    end
    %Calculate interpolative model.
    model = fit(sampled_frequencies(:),calibration_coefficients(:,calibration_idx),'linear');
    save(['speaker_model_' int2str(round(calibration_sampling_frequency/1000))],'model');
    if plot_interpolation
        figure;
        plot(spectrogram_range(1):100:spectrogram_range(2),feval(model,spectrogram_range(1):100:spectrogram_range(2)))
        hold on
        plot(sampled_frequencies,calibration_coefficients(:,calibration_idx),'r.')
        hold off
        xlabel('Frequency (Hz)');
        ylabel('Conversion Coefficient');
        title({'Voltage Power to Pressure Power Conversion Model',...
               ['for ' int2str(round(calibration_sampling_frequency/1000)) 'kHz Sampling Frequency']});
        set(gca,'YScale','log');
        ylim([10^-5 1])
    end
end
