function SetupMenubar(handle)
%edits/deletes items from the standard toolbar and Menus
%handle must be the figure handle

%Author: Amy Boyle
%Last modified June 2010

    %get(findall(handle),'tag')

    %Deletes unnessesary items from the toolbar and menu
    delete([findall(handle, 'tag', 'Standard.NewFigure')...
        findall(handle, 'tag', 'Exploration.Pan') ...
        findall(handle, 'tag', 'Exploration.Rotate') ...
        findall(handle, 'tag', 'Annotation.InsertLegend')...
        findall(handle, 'tag', 'Annotation.InsertColorbar')...
        findall(handle, 'tag', 'DataManager.Linking')...
        findall(handle, 'tag', 'Exploration.Brushing')...
        ])

    delete([findall(handle, 'tag', 'figMenuDesktop')...
        findall(handle, 'tag', 'figMenuHelp') ...
        findall(handle, 'tag', 'figMenuWindow')...
        findall(handle, 'tag', 'figMenuRotate3D')...
        findall(handle, 'tag', 'figMenuPan')...
        findall(handle, 'tag', 'figMenuInsertLegend')...
        findall(handle, 'tag', 'figMenuInsertColorbar')...
        findall(handle, 'tag', 'figMenuCameraToolbar')...
        findall(handle, 'tag', 'figLinked')...
        findall(handle, 'tag', 'figDataManagerBrush')...
        findall(handle, 'tag', 'figBrush')...
        ])
    
%     set(findall(handle, 'tag', 'figMenuFile'), 'visible', 'off');
%     
%     tbh = findall(handle, 'tag', 'FigureToolBar');
%     %get(findall(handle, 'tag', 'figMenuFileExportSetup'), 'callback')
%     
%     set(tbh, 'tag', 'MyToolBar')
%     
%     saveH = findall(handle, 'tag', 'Standard.SaveFigure');
%     set(saveH, 'clickedcallback', @saveProj, 'tooltipstring', 'Save');
%     
%     openH = findall(handle, 'tag', 'Standard.FileOpen');
%     set(openH, 'clickedcallback', 'OpenExisting');
% 
%     fileMenu = uimenu(handle, 'label', 'File', 'position', 1);
%     menuOpen = uimenu(fileMenu, 'label', 'Open', 'callback', 'OpenExisting');
%     menuSave = uimenu(fileMenu, 'label', 'Save', 'callback', @saveProj);
%     menuSaveAs = uimenu(fileMenu, 'label', 'Save As...', 'callback', @saveProjAs);
%     menuExport = uimenu(fileMenu, 'label', 'Export Figure...', 'callback', @exportFun);
%     menuClose = uimenu(fileMenu, 'label', 'Close Figure', 'callback', @closeFig);
%     menuCloseAll = uimenu(fileMenu, 'label', 'Close All', 'callback', @closeAll);
% 
%     
%         function saveProjAs(hObject, eventdata)
%             %saves all the figures currently open
%             [fileName, pathName] = uiputfile;
%             if fileName ~= 0
%                 handArray = findobj('type', 'figure');
%                 hgsave(handArray, [pathName fileName], 'all')
%                 fileNamePath = [pathName fileName];
%                 set(handle, 'userdata', fileNamePath);
%             end
%         end
% 
%         function saveProj(hObject, eventdata)
%             fnp = get(handle, 'userdata');
%             if(isempty(fnp))
%                 saveProjAs
%             else
%                handArray = findobj('type', 'figure');
%                hgsave(handArray, fnp, 'all')
% 
%             end
%         end
% 
%         function closeFig(hObject, eventdata)
%             close gcf;
%         end
% 
%         function closeAll(hObject, eventdata)
%             close all;
%         end
%     
%     function exportFun(hObject, eventdata)
%         %This allows the saving of individual figures
%         filemenufcn(handle,'FileExportSetup')
%     end

end