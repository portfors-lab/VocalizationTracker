function checkBoxSelector
close all
clear variables

import javax.swing.tree.TreeSelectionModel;

%[fileNames pathName] = uigetfile('*.*','Select data files','MultiSelect', 'on')

[mtree, container] = uitree('v0', 'Root','/', 'ExpandFcn', @expandFun); % Parent is ignored
%set(container, 'Parent', hPanel);  % fix the uitree Parent
%mtree.setMultipleSelectionEnabled(true);
jtree = mtree.getTree;

jtree.getSelectionModel().setSelectionMode(TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION)


buttonOne = uicontrol('style', 'pushbutton','position', [255 20 100 40], 'string', 'Load File', 'callback', @loadFun);

%     function nodes = selectFcn(tree, value)
%         nodes = tree.getSelectedNodes;
%         pathList = {};
%         for x = 1:length(nodes)
%             path = nodes(x).getPath.cell;
%             subPathStrs = cellfun(@(p) [p.getName.char, filesep], path, 'un',0);
%             pathStr = strrep([subPathStrs{:}], [filesep,filesep], filesep);
%             pathList = [pathList {pathStr}]
%         end
% 
%     end


     function nodes = expandFun(tree, value)

      try
          count = 0;
          ch = dir(value);
          for i=1:length(ch)
              if (any(strcmp(ch(i).name(1), {'.', ''})) == 0)
                  count = count + 1;
                  if ch(i).isdir
                      iconpath = [matlabroot, '/toolbox/matlab/icons/foldericon.gif'];
                  else
                      iconpath = [matlabroot, '/toolbox/matlab/icons/pageicon.gif'];
                  end
                  nodes(count) = uitreenode('v0',[value, ch(i).name, filesep], ...
                      ch(i).name, iconpath, ~ch(i).isdir);
              end
          end
      catch
          error('MyApplication:UnrecognizedNode', ...
            ['The uitree node type is not recognized. You may need to ', ...
            'define an ExpandFcn for the nodes.']);
      end

      if (count == 0)
          nodes = [];

      end
    end

    function loadFun(hobject, eventdata)
        nodes = mtree.getSelectedNodes;
        pathList = {};
        for x = 1:length(nodes)
            path = nodes(x).getPath.cell;
            subPathStrs = cellfun(@(p) [p.getName.char, filesep], path, 'un',0);
            pathStr = strrep([subPathStrs{:}], [filesep,filesep], filesep);
            pathList = [pathList {pathStr}]
        end
        
        %Create a preferences structure for the desired experimental data
        prefs = GeneratePreferencesNew(pathList{1})
        %Set the threshold used for spike detection. 0.11 is the default.
        prefs.spike_time_peak_threshold = 0.2;
        %Extract XML metadata and conver to to Matlab structure
        experiment_data = LoadExperimentData(prefs);

        %-------------------
        % Test Visualization
        %-------------------
        %Specify the test number to visualize, this is a one tone test
        test_num = 9;
        %Generate contour plot of single frequency tuning curve
        VisualizeTestData(experiment_data,prefs,test_num);
        %Generate image map plots of time-frequency histograms
        VisualizeTestData(experiment_data,prefs,test_num,[0 1 0 0 1]);

        % --------------------
        % Trace Visualization
        % --------------------
        % Specify the test and trace number to visualize
        test_num = 9;
        trace_num = 6; 
        %Visualize the traces
        VisualizeTraceData(experiment_data,prefs,test_num,trace_num);
        VisualizeTraceData(experiment_data,prefs,test_num,trace_num,[0 1 1 0 0 0 0 0 0 0]);
        VisualizeTraceData(experiment_data,prefs,test_num,trace_num,[0 0 0 1 1 1 1 0 0 0]);
        VisualizeTraceData(experiment_data,prefs,test_num,trace_num,[0 0 0 0 0 0 0 1 1 0]);
    end
end
% selNodes = tree.selectedNodes;
% selNode = selNodes(1);
% selNodeName = char ( selNode.getName() )
% path = selNode.getPath;
% parentNode = char(path( (length(path))-1 ).getName) % getting the
% parent of the selected node 