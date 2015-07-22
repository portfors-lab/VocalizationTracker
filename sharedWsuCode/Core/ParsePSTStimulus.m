function [stimulus stim_type]= ParsePSTStimulus(stim_text,test_num,trace_num)

%Get the global stimulus parameters
[stim_data position] = textscan(stim_text,'%d %d %f %f %f',1);

%Check to see if the channel is activated
if stim_data{1} == 0
    stimulus = [];
    stim_type = 0;
    return;
end
%Set the global stimulus parameters
stim_type = stim_data{2};
stimulus.attenuation = stim_data{3};
stimulus.duration = stim_data{4};
stimulus.delay = stim_data{5};

%Chop off the part of the stimulus string already read above
stim_text = stim_text(position+1:end);

%Static stimulus labels
stim_types = {'tone',...
              'fmsweep',....
              'synthesized_batsound',....
              'amsound',....
              'broad_band_noise',....
              'narrow_band_noise',....
              'click',....
              'stored_vocal_call',....
              'high_pass_noise',....
              'low_pass_noise',....
              'sine_wave_modulation',....
              'square_wave_modulation'};

%Default values
stimulus.frequency = [];
stimulus.rise_fall = [];
stimulus.soundtype_name = [];
stimulus.reverse_vocal_call = [];
stimulus.vocal_call_file = [];
stimulus.bandwidth=[];
stimulus.usweep=[];

switch stim_type
    case 0,
        return
    case 1, %Tone
        stim_data = textscan(stim_text,'%f %f',1);
        stimulus.frequency = stim_data{1};
        stimulus.rise_fall = stim_data{2};
        stimulus.soundtype_name = 'tone';
    case 2
        stim_data = textscan(stim_text,'%f %f %f %f',1);
        stimulus.frequency=stim_data{1};
        stimulus.bandwidth=stim_data{2};
        stimulus.usweep=stim_data{3};
        stimulus.rise_fall = stim_data{4};
        stimulus.soundtype_name = 'fmsweep';
    case 8, %Vocaization
        stim_data = textscan(stim_text,'%d',2);
        stim_data = stim_data{1};
        %Note: Lots of variables ignored. This would be a good place to add
        %stuff if it is found that modified vocalizations are used and that
        %the modifications need to be accounted for in Bat2Matlab
        stimulus.reverse_vocal_call = stim_data(1);
        stimulus.soundtype_name = 'vocalization';
        if stim_data(2) == 0
            %Old school
            stim_data = textscan(stim_text,'%s',16);
            stim_data = stim_data{1};
            stimulus.vocal_call_file = stim_data{16};
        else
            %New School
            stim_data = textscan(stim_text,'%s',20);
            stim_data = stim_data{1};
            if length(stim_data{20}) > 1 && ~strcmp(stim_data{20},'-1')
                %*PDR* error('PST file generated between March and September 2007. No audio file extension information available.');
                warning('PDR:changes', 'PST file generated between March and September 2007. \n No audio file extension information available. \n Using .call1');
                stimulus.vocal_call_file = [stim_data{20} '.call1']; %*PDR*12/04/07
            else
            %Get vocalization audio file information
            %First, skip ahead to the vocalization info
                [stim_data position] = textscan(stim_text,'%d',28);
                stim_text = stim_text(position+1:end);
                [stim_data position] = textscan(stim_text,'%s',6);
                stim_data = stim_data{1};
                stimulus.vocal_call_file = [stim_data{1} stim_data{3}];
            end
        end
    otherwise,
        stimulus.soundtype_name = stim_types{stim_type};
%         warning(['Unsupported Stimulus Type ' stim_types{stim_type} '. Stimulus not imported for test ' int2str(test_num) ', trace ' int2str(trace_num) '.']);
end
