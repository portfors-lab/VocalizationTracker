kHz = 10^3;
% Specify the audio file
path = '';
fname = 'Am_sm6.call1';
filePath = [path '\' fname];
% Load the audio file
[y sampleRate] = ParseAudioData(filePath);
%Filter out hum at 70 kHz
notchCenter = 70*kHz/(sampleRate/2);
notchWidth = notchCenter/300;
[b,a] = iirnotch(notchCenter,notchWidth);  
y = filtfilt(b,a,y);
