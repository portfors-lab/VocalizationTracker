function [signal sampleRate interval] = ParsePartialAudio(filePath)
%get portion of the wav file returning the signal, sample rate and
%interval time in seconds

%Amy Boyle 6/2/11
   
    maxNSamples = 1e7; %fairly arbitrary, change as fits
    
    if strfind(filePath, 'wav')
        [signal sampleRate] = wavread(filePath, [1 10]); %just to get sample rate
        [ws wavstring] = wavfinfo(filePath);
        if isempty(ws)
            warndlg('only longer wav files permitted');
            signal = [];
            return
        end
        nsamples = textscan(wavstring, 'Sound (WAV) file containing: %u');
        nsamples = nsamples{1};
        seconds = nsamples/sampleRate; %duration of file
        maxTime = maxNSamples/sampleRate; %longest chunk of time possible to view
        interval = timeRangeDlg;
        if isempty(interval)
            signal = [];
            return
        end
        while (interval(2)-interval(1))*sampleRate > maxNSamples
            interval = timeRangeDlg('interval exceeds limit, try again');
            if isempty(interval)
                signal = [];
                return
            end
        end
        if interval(1) ==0
            start = 1;
        else
            start = interval(1)*sampleRate;
        end
        stop = interval(2)*sampleRate;
        signal = wavread(filePath, [start stop]);
        signal = signal(:,1); %if stereo, get only left channel
        signal = signal./max(abs(signal)); %normalize
    elseif strfind(filePath, 'call') || strfind(filePath, 'kanwal')
        if strfind(filePath, 'call')
            sampleRate = 333333; %always assumed for call files
        else
            sampleRate = 250000; %always assumed for kanwal files
        end
        finfo = dir(filePath);
        nsamples = finfo.bytes/2;
        seconds = nsamples/sampleRate; %duration of file
        maxTime = maxNSamples/sampleRate; %longest chunk of time possible to view
        interval = timeRangeDlg;
        if isempty(interval)
            signal = [];
            return
        end
        while (interval(2)-interval(1))*sampleRate > maxNSamples
            interval = timeRangeDlg('interval exceeds limit, try again');
            if isempty(interval)
                signal = [];
                return
            end
        end
        start = interval(1)*sampleRate; %in samples
        stop = interval(2)*sampleRate;
        duration = stop-start;
        fid = fopen(filePath, 'r', 'l');
        fseek(fid,start*2,-1); %fseek uses bytes, so we multiply x2 as 2 bytes/sample
        signal = fread(fid,duration,'int16'); % we define precision to be int16, so duration is in samples.
    else
        disp('Unsupported file type for longer duration files');
    end
    
    function range = timeRangeDlg(errormsg)
        range = [];
        fh = figure('position', [400 250 500 225], 'name', 'Select time interval', 'windowstyle', 'modal', 'numbertitle', 'off', 'menubar', 'none');
        msg = ['Select time interval between 0 and ' num2str(seconds) 'seconds, not to exceed ' num2str(maxTime) ' seconds'];
        uicontrol(fh, 'style', 'text', 'position', [50 160 400 30], 'string', msg, 'backgroundcolor', [0.8 0.8 0.8]);
        startField = uicontrol(fh, 'style', 'edit', 'position', [100 100 100 30]);
        stopField = uicontrol(fh, 'style', 'edit', 'position', [300 100 100 30]);
        uicontrol(fh, 'style', 'pushbutton', 'position', [300 25 75 30], 'string', 'ok', 'callback', @okFun)
        uicontrol(fh, 'style', 'pushbutton', 'position', [400 25 75 30], 'string', 'cancel', 'callback', 'close(gcf)')
        if exist('errormsg', 'var')
            uicontrol(fh, 'style', 'text', 'position', [50 150 400 20], 'string', errormsg, 'foregroundcolor', 'red', 'backgroundcolor', [0.8 0.8 0.8]);
        end
        uiwait
        
        function okFun(ho,ed)
            range(1) = str2double(get(startField, 'string'));
            range(2) = str2double(get(stopField, 'string'));
            if any(isnan(range))
                range = [];
            end
            close(fh)
        end
    end
end