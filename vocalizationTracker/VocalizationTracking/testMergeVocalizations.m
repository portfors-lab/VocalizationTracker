close all
clear variables

path = '.';
fname = 'Am_cm12.call1';
filePath = [path '/' fname];
[y sampleRate] = ParseAudioData(filePath);

% Put all of the different vocalization components into a cell array
vocalizationCellArray = {y,y,y};
mergedVocalization = MergeVocalizations(vocalizationCellArray, sampleRate);
figure;
plot(mergedVocalization);
