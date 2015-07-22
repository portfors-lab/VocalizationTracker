function [signal sample_rate] = ParseAudioData(file_path)
%
%[signal sample_rate] = ParseAudioData(file_path)
%
%   INPUT ARGUMENTS
%   file_path       Path to the audio file.
%
%   OUTPUT ARGUMENTS
%   signal          Audio signal, normalized to 1.
%   sample_rate     The sample rate of the audio signal

if strfind(file_path,'.kanwal')
    fid = fopen(file_path,'r','l'); %Open the file
    if fid == -1
        error(['Could not find vocalization file located at path: ' file_path]);
    end
    signal = fread(fid,inf,'int16')';
    fclose(fid);
    sample_rate = 2.5e5; %250kHz sample rate
elseif strfind(file_path,'.call')
    fid = fopen(file_path,'r','l'); %Open the file
    if fid == -1
        error(['Could not find vocalization file located at path: ' file_path]);
    end
    signal = fread(fid,inf,'int16')';
    fclose(fid);
    sample_rate = 333333; %333,333Hz sample rate
elseif strfind(file_path,'.wav')
    fid = fopen(file_path,'r','l'); %Open the file
    if fid == -1
        error(['Could not find vocalization file located at path: ' file_path]);
    end
    %Rewind file pointer
    fseek(fid,0,-1);
    %Grab first chunk of file to look for start of wav file
    S = fscanf(fid,'%c',3000);
    RIFF_offsets = strfind(S,'RIFF');
    if length(RIFF_offsets) == 1
        %Standard wav file format
        fclose(fid);
        [signal sample_rate] = wavread(file_path);
        %Grab left channel if it is a stereo recording
        signal = signal(:,1);
    elseif length(RIFF_offsets) == 3
        %Extract sample rate out of file
        fseek(fid,hex2dec('104'),-1);
        sample_rate = fread(fid,1,'uint32') * 2;
        %Location of 3rd 'RIFF' location plus 44 byte wav header minus one for
        %zero based hex indexing
        wav_start = RIFF_offsets(3)+43;
        %Scan to start of wav file
        fseek(fid,wav_start,-1);
        signal = fread(fid,inf,'int16')';
        fclose(fid);
    else
        error('wav file is corrupted.');
    end

else
    error('Unsupported Audio File Type.');
end

% %Normalize
signal = signal./max(abs(signal));
%% --- High-pass filter ---
% Hd = highpass10khz;
% % Hd = highpass30khz;
% signal = filter(Hd,signal);
% %% --- [Lukashkin98] cochlea ---
% scale = 20;
% signal = Lukashkin(signal*scale, sample_rate);
%% ------------------------

