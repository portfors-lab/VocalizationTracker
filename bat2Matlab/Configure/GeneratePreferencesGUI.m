function prefs = GeneratePreferencesGUI(animal_type, ...
                                     animal_string, ...
                                     experiment_string, ...
                                     cell_string)
%
%function prefs = GeneratePreferences(animal_number,
%                                     experiment_letter,
%                                     cell_depth)
%
%   INPUT ARGUMENTS
%   animal_type         The type of animal used in the experiment.
%                       Example: 'bat' or 'mouse'
%   animal_string       The descriptor of the animal used for the experimental
%                       data, in string format.
%                       Example: '7'
%   experiment_string   The descriptor of the specific experiment used for the
%                       experimental data, in string format.
%                       Example: 'c'
%   cell_string         The desciptor of the cell used for the experimental data. 
%                       The electrode depth is a good choice.
%                       Example: '1440'
%
%   OUTPUT ARGUMENTS
%   prefs               The preferences used throughout the Bat2Matlab
%                       execution.
%
%GeneratePreferences generates a structure containing all of the
%global preferences used throughout the execution of Bat2Matlab.

if exist('animal_type','var')
    %Test description
    prefs.animal_type = animal_type;
    prefs.animal_string = animal_string;
    prefs.experiment_string = experiment_string;
    prefs.cell_string = cell_string;

    %Path definitions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    base_batlab_data_path = uigetdir('C:\LAB RESEARCH\Bat2Matlab', 'Select a directory');
prefs.bat2matlab_directory = base_batlab_data_path;
%     if strmatch('bat',animal_type)
%         prefs.audio_directory = [prefs.bat2matlab_directory '\Data\'];
%         base_batlab_data_path = [prefs.bat2matlab_directory '\Data\'];
%     elseif strmatch('mouse',animal_type)
%         prefs.audio_directory = [prefs.bat2matlab_directory '\Data'];
%         base_batlab_data_path = [prefs.bat2matlab_directory '\Data'];
%     else
%         error('Incorrect animal type specified in GeneratePreferences()');
%     end

 [tok rest] = strtok(base_batlab_data_path, '\');
parentDir = [tok];
 while ~isempty(rest)
 [tok rest] = strtok(rest, '\')
 parentDir = [parentDir '\' tok]
end
parentDir = parentDir(1:length(parentDir) - length(tok))

    prefs.Bat2Matlab_data_filepath = [parentDir '\Extracted Data\' tok '.mat'];
    prefs.raw_data_filepath = [base_batlab_data_path '\' tok '.raw'];
    prefs.xml_data_filepath = [base_batlab_data_path '\' tok '-alltests-nospikes.xml'];
    prefs.pst_data_filepath = [base_batlab_data_path '\' tok '.pst'];
    prefs.output_data_filepath = [prefs.bat2matlab_directory '\Output\'];
    if ~isempty(cell_string)
        prefs.output_data_filepath = [prefs.output_data_filepath '_' cell_string];
    end
    prefs.cache_dir = [prefs.output_data_filepath '\cache'];
    prefs.cell_id = [tok '_' cell_string];
    prefs.cell_id4_plot = prefs.cell_id; prefs.cell_id4_plot(strfind(prefs.cell_id4_plot,'_')) = '.';
end

prefs.speaker_calibration_file = 'speaker_calibration_8_27_06';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options used for speaker callibration and dB SPL normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20 micro-Pascals for dB SPL conversion
prefs.dB_SPL_ref = 0.00002;
%16-bit audio dynamic range
prefs.dbRange = 96.3296;
%Approximate maximum output of speaker in dB SPL
prefs.dbMax = 115;
%Approximate maximum output of speaker in dB SPL minus maximum dynamic range of 16 bit audio.
% prefs.dbMin = dbMax-dbRange;
prefs.dbMin = 0;
%This allows the setting of a threshold below which everything is set to zero.
prefs.model_data_dbMin = prefs.dbMin;
%Approximate maximum output of speaker in dB SPL/Hz
prefs.spectrogram_dbMax = 84.5665;
%Approximate maximum output of speaker in dB SPL/Hz minus maximum dynamic range of 16 bit audio.
prefs.spectrogram_dbMin = prefs.spectrogram_dbMax-prefs.dbRange;
% prefs.spectrogram_dbMin = 0;
%Decide whether to use absolute or relative color scaling when plotting the spectrogram.
prefs.spectrogram_absolute_scaling = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Spectrogram calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The number of times to evaluate the spectrogram in the frequency domain.
prefs.spectrogram_freq_samples = 2^10;
%The density of sampling in the time domain
% prefs.spectrogram_time_samples_per_millisecond = 1/2;
% prefs.spectrogram_time_samples_per_millisecond = 2;
prefs.spectrogram_time_samples_per_millisecond = 5;
%Frequency range for spectrogram
prefs.spectrogram_range = [0 110000];
%Window length (in seconds)
% prefs.spectrogram_window_length = 0.0001;
% prefs.spectrogram_window_length = 0.0005;
% prefs.spectrogram_window_length = 0.0015; %Happy medium? Default
prefs.spectrogram_window_length = 0.002; %Best so far
% prefs.spectrogram_window_length = 0.0025;
% prefs.spectrogram_window_length = 0.003;
% prefs.spectrogram_window_length = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Spike Time calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs.spike_time_filter_cutoff = 600;
prefs.spike_time_power_exponent = 2;
prefs.spike_time_peak_threshold = 0.11; %Default
prefs.spike_time_refractory_period = 2.0; %Milliseconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Spike Rate calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prefs.spike_rate_sampling_frequency = NaN; %defaults to trace.samplerate_ad;
prefs.spike_rate_sampling_frequency = 8000;
prefs.spike_rate_cutoff_frequency = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for spectrogram peak finding calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs.num_harmonics = 6;
%Peak detection noise floor as a fraction of the highest peak
prefs.harmonic_peak_detection_noise_floor = 0.2;
prefs.frequency_cutoff_ratio = 175; %Higher values smooth more

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for model data generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs.model_num_data_rows_per_cell = 2^7;
prefs.model_time_samples_per_millisecond = 1/2; %Default
% prefs.model_time_samples_per_millisecond = 1; %Default
prefs.model_spectral_integration = 0; %0:Rectangular 1:Gaussian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Histogram generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prefs.histogram_bin_width = 1; 
prefs.histogram_bin_width = 4; %Default
% prefs.histogram_bin_width = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for the calculation of the frequency intervals to integrate over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs.max_interval_width = 4000;
prefs.default_intervals = [10000:2000:98000 ; 12000:2000:100000]';
prefs.default_sampled_frequencies = [11000:2000:99000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for harmonic peak calculation as a fraction of the highest peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for VR metric calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs.exponential_decay = 10; %Milliseconds
%The sample frequency of the signal reconstituted from thespike times
prefs.filtered_exponential_sample_frequency = 1000; %Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Spike Train Filtering (Gaussian kernal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prefs.filtered_gaussian_sample_frequency = 1000; %Hz
% prefs.filtered_gaussian_stdev = 2; %In milliseconds
prefs.filtered_gaussian_stdev = 3; %In milliseconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global options for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Calling the colormap functions creates a figure if there is not
%already one created. Here we check to see if a figure exists already.
%If not, delete the figure generated by the colormap function.
figure_exists = get(0,'CurrentFigure');
prefs.force_plot_visible = false;
prefs.colormap = jet;prefs.colormap_name = 'jet';
% prefs.colormap = hot;prefs.colormap_name = 'hot';
% prefs.colormap = gray;prefs.colormap_name = 'gray';
% prefs.colormap = gray;prefs.colormap_name = 'bone';
% prefs.colormap = colorspiral;prefs.colormap_name = 'colorspiral';
% prefs.colormap = flipud(prefs.colormap); %Use this to reverse colormap
% prefs.colormap = brighten(prefs.colormap,0.3);
if isempty(figure_exists)
    close;
end

    