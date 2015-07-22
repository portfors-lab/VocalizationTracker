function [prefs experimentData] = GetExpData(animalPath, isAnimalName)

%checks to see if experimentData is cached in a .mat file before parsing
%.pst and .raw files.

    if exist('isAnimalName', 'var')
        %given the animal name, find files that end with .pst and .raw in the
        %same folder.
        [~, animalName] = fileparts(animalPath);
        animalMatFile = [animalName '_cachedData.mat'];
        prefs=GeneratePreferencesNew2(animalPath);
        cachedData=[prefs.extracted_data_folder animalMatFile];
        if exist(cachedData, 'file')
            d = load(cachedData, '-mat');
            experimentData=d.experimentData;
        else
            if ~exist(prefs.raw_data_filepath, 'file') || ~exist(prefs.pst_data_filepath, 'file') 
                WriteStatus('Invalid animal name, select the name of the files which end in .pst and .raw (but leave extention off)', 'red');
                prefs = [];
                experimentData = [];
                return
            end 
            experimentData = LoadExperimentData(prefs);
            save(cachedData, 'experimentData');
        end

    else %older way of doing it where the folder name and file names match
        animalMatFile = 'cachedData.mat';
        if ~isequal(animalPath(end), filesep)
            animalPath = [animalPath filesep];
        end
        cachedData = [animalPath animalMatFile];
        prefs = GeneratePreferencesNew(animalPath);
        if exist(cachedData, 'file')
            d = load(cachedData, '-mat');
            experimentData = d.experimentData;
        else
            if ~exist(prefs.raw_data_filepath, 'file') || ~exist(prefs.pst_data_filepath, 'file')
                WriteStatus('Invalid animal folder selection: select the folder in which the .raw & .pst files reside; folder must bear the same name as these files', 'red');
                prefs = [];
                experimentData = [];
                return
            end
            experimentData = LoadExperimentData(prefs);
            save(cachedData, 'experimentData');
        end
    end
end