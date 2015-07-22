function WriteAudioData(signal,sample_rate,file_path)
%
%WriteAudioData(signal,sample_rate,file_path)
%
%   INPUT ARGUMENTS
%   signal          The signal to be converted
%   sample_rate     The sample rate of the signal to be converted
%   file_path       Path to the target audio file, including the extension.
%                   File extensions include 'kanwal', 'call', , 'call1', and 'wav'

%Normalize
signal = signal./max(signal);

if strfind(file_path,'.kanwal')
    if sample_rate ~= 2.5e5
        %resample cannot handle the large sample rate numbers so these are
        %equivalent ratios for common sampling rates.
        warning('Resampling, file quality being reduced');
        switch(sample_rate)
            case(500000)
                signal = downsample(signal,2);
            case(333333)
                signal = resample(signal, 3,4);
            otherwise
                error('Sampling rate non-standard')
        end
    end
    signal = PCM_Quantize(signal);
    fid = fopen(file_path,'w'); %Create the file
    fwrite(fid,signal,'int16');
    fclose(fid);
elseif strfind(file_path,'.call')
    if sample_rate ~= 333333
        %error('Call audio files need to be 333,333 Hz');
        warning('Resampling, file quality being reduced');
        switch(sample_rate)
            case(500000)
                signal = resample(signal,2,3);
            case(250000)
                signal = resample(signal,4,3);
            case(375000)
                signal = resample(signal,8,9);
            otherwise
                error('Sampling rate non-standard')
        end
    end
    signal = PCM_Quantize(signal);
    fid = fopen(file_path,'w'); %Create the file
    fwrite(fid,signal,'int16');
    fclose(fid);
elseif strfind(file_path,'.wav')
    warning off
    wavwrite(signal,sample_rate,file_path)
    warning on
else
    error('Unsupported Audio File Type.');
end

% -----------------------------------------------------------------------
function y = PCM_Quantize(x)
% PCM_Quantize:
%   Scale and quantize input data, from [-1, +1] range to 16-bit data range

% Determine slope (m) and bias (b) for data scaling:
nbits = 16;
m = 2.^(nbits-1);
b = 0;

y = round(m .* x + b);

% Determine quantized data limits, based on the
% presumed input data limits of [-1, +1]:
ylim = [-1 +1];
qlim = m * ylim + b;
qlim(2) = qlim(2)-1;

% Clip data to quantizer limits:
i = find(y < qlim(1));
if ~isempty(i),
   y(i) = qlim(1);
end

i = find(y > qlim(2));
if ~isempty(i),
   y(i) = qlim(2);
end

return