spectrogram_range = [0 100000];
%The results of the callibration test
sampled_frequencies = 6:2:100;
max_dB_SPL = [99 94 100.5 108.5 111.0 108.5 104.5 102 100.5 98 97.5 99.5 100.5 101.5 98.5 96 93 92 91 92.5 92.5 91 90.5 89.5 88 87 86 85 84.5 83.5 83 82 83.5 82.5 82.5 83 83 80 80 79 79 79 79 79 79.5 79 79.5 80];
%Add samples before and after to constrain extrapolation
pre_sampled_frequencies = [2 4];
pre_max_dB_SPL = ones(1,length(pre_sampled_frequencies))*max_dB_SPL(1);
post_sampled_frequencies = [];
post_max_dB_SPL = ones(1,length(post_sampled_frequencies))*max_dB_SPL(end);
sampled_frequencies = [pre_sampled_frequencies sampled_frequencies post_sampled_frequencies]*1000;
max_dB_SPL = [pre_max_dB_SPL max_dB_SPL post_max_dB_SPL];

%For testing. Calculate interpolative model and plot interpolations.
[curve, goodness] = fit(sampled_frequencies(:),max_dB_SPL(:),'linear');
figure;
plot(spectrogram_range(1):100:spectrogram_range(2),feval(curve,spectrogram_range(1):100:spectrogram_range(2)));
hold on
plot(sampled_frequencies,max_dB_SPL,'r.')
hold off
xlabel('Frequency (Hz)');
ylabel('Speaker Output in dB SPL');
title('Speaker Calibration Curve');