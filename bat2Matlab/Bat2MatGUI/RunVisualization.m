function RunVisualization(javaPanel, varargin)

    errorCheckFirst = false; %whether to check for user input error initally or as you go

    saveOn = false;
    txtPath = [];
    warnColor = 'systemcommands';
    err = 'red';
    warn = 'orange';
    fig = gcf;
    info = 'black';
    savePath =[];
    
    if exist('varargin', 'var')
        if mod(length(varargin),2) ~= 0
            error(['Unexpected number of arguments. Optional arguments ' ...
                      'must be specified as property value paris.']);
        end
        for ind = 1:2:length(varargin)
            parameterName = varargin{ind};
            switch lower(parameterName)
                case 'saveon'
                    savePath = varargin{ind+1}{1};
                    format = varargin{ind+1}{2};
                    resolution = varargin{ind+1}{3};
                    saveOn = true;
                case 'rootpath'
                    rootPath = varargin{ind+1};
                case 'colormap'
                    colormap =  varargin{ind+1};
                case 'invertcolor'
                    invertColor = varargin{ind+1};
                case 'stimpath'
                    stimPath = varargin{ind+1};
                case 'txtpath'
                    txtPath = varargin{ind+1};
                case 'colorrange'
                    colorRange = varargin{ind+1};
            end

        end
    end

    try
        paths = cellstr(char(javaPanel.getCheckedPaths()));
        if isempty(paths{1})
            WriteStatus('Select a file to do sciency stuff with', err)
            return
        end
        %get what plots to make
        testChoices = javaPanel.getTestChoices();
        traceFlags = javaPanel.getTraceChoices();
        if ~any([testChoices; traceFlags])
%             disp('Silly Fool, make some plot choices!')
            WriteStatus('You must make some plot choices, Ya Dafty', err)
            return
        end
        if testChoices(1) || testChoices (2)
            imFlags = [testChoices(1) testChoices(2) 0 0 1];
        else
            imFlags = 0;
        end
        if testChoices(3) || testChoices(4)
            conFlags = [testChoices(3) testChoices(4) 0 1 0];
        else 
            conFlags = 0;
        end
        if testChoices(5) || testChoices(6)
            surfFlags = [testChoices(5) testChoices(6) 1 0 0];
        else 
            surfFlags = 0;
        end
        
        if errorCheckFirst
            mods = javaPanel.getAllMods();
            while(mods.hasMoreElements());
                module = mods.nextElement();
                traceTests = char(module.getTests2());
                if ~isempty(traceTests)
                    traces = module.getTraces();
                    if isempty(traces)
                        error(['You must select an option for trace numbers for ' char(module.getName())])
                    end
                else
                    tests = char(module.getTests1());
                    if isempty(tests)
                        error(['You must select at least one test/trace for ' char(module.getName())])
                    end
                end
                threshold = str2double(module.getThreshold());
                binSize = str2double(module.getBinsize());
                if isnan(threshold) || isnan(binSize)
                    error(['invalid threshold/bin size for ' char(module.getName())]);
                end                
            end
        end
        
       
        for a = 1:length(paths)
            
            %--------------------------------------------------------------
            % Get data for Experiment run
            %--------------------------------------------------------------
            
            path = [rootPath paths{a}];
            %Create a preferences structure for the desired experimental data
            [prefs experiment_data] = GetExpData(path);
            if ~exist(prefs.extracted_data_folder, 'file')
                mkdir(prefs.extracted_data_folder);
            end
            name = prefs.cell_id(1:end-1);
            WriteStatus(['Processing ' name], info, fig);
            %user defined preferences
            threshold = str2double(javaPanel.getThreshold(name));
            binSize = str2double(javaPanel.getBinsize(name));
            if isnan(threshold)
               WriteStatus(['Human induced malfunction- Invalid theshold value, skipping experiment.'], warn, fig);
               continue
            end
            if isnan(binSize)
                WriteStatus(['Human error! Invalid bin size value, skipping experiment. Silly Human'], warn, fig);
                continue
            end
            prefs.spike_time_peak_threshold = threshold;
            prefs.histogram_bin_width = binSize;
            prefs.audio_directory = stimPath;
            prefs.colormap = eval(colormap);
            prefs.colormap_name = colormap;
            prefs.invert_color = invertColor;
                        
            %==============================================================
            %Test Visualizsation
            %--------------------------------------------------------------
            %get the test numbers to do whole test analysis on
            testNumString = char(javaPanel.getTests1(name));
            if ~isempty(testNumString)
                if ~isempty(txtPath)
                    WriteStatus('Tuning Curves not exportable as text format', warn, fig);
                else
                %get the test numbers and check input for errors
                testNumList = evalStrList(testNumString);
                %produce desired plots for each test
                for testNum = testNumList 
%                     disp(['Processing Test ' num2str(testNum) ' tuning curve...']);
                    WriteStatus(['Processing Test ' num2str(testNum) ' tuning curve...'], info, fig);
                    try
                        if any(imFlags)
                            if(saveOn)
                                VisualizeTestData(experiment_data,prefs,testNum,imFlags, savePath, format, resolution, colorRange);
                            else
                                VisualizeTestData(experiment_data,prefs,testNum,imFlags, [], [], [], colorRange);
                            end                            
                        end
                        
                        if any(conFlags)
                           if(saveOn)
                              VisualizeTestData(experiment_data,prefs,testNum,conFlags, savePath, format, resolution, colorRange);
                           else
                              VisualizeTestData(experiment_data,prefs,testNum,conFlags, [], [], [], colorRange);    
                           end
                        end
                        
                        if any(surfFlags)
                            if(saveOn)
                                VisualizeTestData(experiment_data,prefs,testNum,surfFlags, savePath, format, resolution, colorRange);
                            else
                                VisualizeTestData(experiment_data,prefs,testNum,surfFlags, [], [], [], colorRange);   
                            end
                        end                        
                    catch e
                        WriteStatus('unexpected error occured, see command line', err, fig);
                        cprintf('err', ['Test ' num2str(testNum) ': ' e.message '\n']);
                        %rethrow(e)
                        continue
                    end
                    
                end
                end
            end
            %==============================================================
            %Trace Visualization
            %--------------------------------------------------------------
            %Do trace analysis
            clear testNumString;
            clear traceNumString;
            
            %get test numbers for analysis, do nothing if empty
            testNumString = char(javaPanel.getTests2(name));
            if ~isempty(testNumString)
                testNumList = evalStrList(testNumString);
                
                selectionType =  char(javaPanel.getSelectionType(name));
                traceNumString = char(javaPanel.getTraces(name));
                                
                if isempty(traceNumString)
                    WriteStatus(['You must select a choice for traces'], warn,fig);
                %all radio button selected 
                elseif strcmp(selectionType, 'all')
                    if ~isempty(txtPath)
                        TextFileOutput(experiment_data, prefs, testNumList, [], traceFlags, txtPath);
                        continue
                    end
                    for testNum = testNumList
%                         disp(['Processing Test ' num2str(testNum) ' trace data...']);
                        WriteStatus(['Processing Test ' num2str(testNum) ' trace data...'], info, fig);
                        %get all the traces for each test
                        traceNums = 1:length(experiment_data.test(testNum).trace);
                        try
                            for traceNum = traceNums                               
                                if saveOn
                                    VisualizeTraceData(experiment_data,prefs,testNum,traceNum,traceFlags, savePath, format, resolution);
                                else
                                    VisualizeTraceData(experiment_data,prefs,testNum,traceNum,traceFlags);   
                                end
                            end
                        catch e
                            WriteStatus('unexpected error occured, see command line. Continuing...', err,fig);
                            cprintf('err', ['Test ' num2str(testNum) ': ' e.message '\n']);
                            continue
                            %rethrow(e)
                        end
                    end
                elseif strcmp(selectionType, 'same')
                    traceNums = str2num(traceNumString);
                    if isempty(traceNums)
                        WriteStatus('Invalid trace number input', warn,fig);
                        continue
                    end
                    if ~isempty(txtPath)
                        TextFileOutput(experiment_data, prefs, testNumList, traceNums, traceFlags, txtPath);
                        continue
                    end
                     for testNum = testNumList
                        WriteStatus(['Processing Test ' num2str(testNum) ' trace data...'], info, fig);
                        for thisTrace = traceNums 
                            try
                                if saveOn
                                    VisualizeTraceData(experiment_data, prefs, testNum, thisTrace,traceFlags,savePath, format,resolution);
                                else
                                    VisualizeTraceData(experiment_data, prefs, testNum, thisTrace,traceFlags);  
                                end
                            catch e
                                WriteStatus('unexpected error occured, see command line. Continuing...', err,fig);
                                cprintf('err', ['Test ' num2str(testNum) ': ' e.message '\n']);
                                continue
                                %rethrow(e)
                            end
                        end
                     end
                %select traces radio selected
                elseif strcmp(selectionType, 'different')
                    
                    traceCell = textscan(traceNumString, '%s', 'delimiter', ';');
                    traceCell = cellfun(@str2num, traceCell{1}, 'uniformoutput', false);
                    traceIdx = 1;
                    if ~isempty(txtPath)
%                         traceNumArray  = evalStrList2(traceNumString);
                        TextFileOutput(experiment_data, prefs, testNumList, traceCell, traceFlags, txtPath);
                        continue
                    end
                    
                    for testNum = testNumList
                        WriteStatus(['Processing Test ' num2str(testNum) ' trace data...'], info, fig);
                        %matches each test num to its trace num(s) in the
                        %trace num field
                        traceNum = traceCell{traceIdx};
                        traceIdx = traceIdx +1;
                        if isempty(traceNum)
                            WriteStatus('Invalid trace number input',warn,fig);
                            continue
                        end
                        for thisTrace = traceNum                                        
                            try
                                if saveOn
                                    VisualizeTraceData(experiment_data, prefs, testNum, thisTrace,traceFlags,savePath, format,resolution);
                                else
                                    VisualizeTraceData(experiment_data, prefs, testNum, thisTrace,traceFlags);  
                                end
                            catch e
                                WriteStatus('unexpected error occured, see command line. Continuing...', err,fig);
                                cprintf('err', ['Test ' num2str(testNum) ': ' e.message '\n']);
                                continue
                                %rethrow(e)
                            end
                        end
                    end
                end
            end
        end
    catch e
    %             set(hobject, 'enable', 'on');
%             errordlg(e.message)
            cprintf('err', [e.message '\n']);
            WriteStatus('unexpected error, check command line', err,fig);
            rethrow(e)
    end
end