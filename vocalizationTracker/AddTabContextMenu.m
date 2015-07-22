function AddTabContextMenu(thisTab, varargin)
%adds a context menu to plots residing in tabs

%Amy Boyle  3/24/11
%          12/28/11

    if ~isempty(varargin)
        if rem(length(varargin),2)~=0
           disp('optional Arguments to AddTabContextMenu must be in name value pairs');
           return
        end
        for indx=1:2:length(varargin)-1
            switch(lower(varargin{indx}))
                case 'saveplot'
                    saveplot=varargin{indx+1};
                case 'delete'
                    deleteplot=varargin{indx+1};
                case 'changecolormap'
                    changecmap=varargin{indx+1};
            end
        end
    else
        saveplot=true;
        deleteplot=true;
        changecmap=true;
    end

    %Must create a context menu object for each tab, if you only create
    %one and assign it to all, it gets deleted with the tab if you
    %delete the tab
    hTabGroup = get(thisTab, 'parent');
    cmenu = uicontextmenu;
    if saveplot
        uimenu(cmenu, 'label', 'save plot', 'callback', @saveTabbedPlot);
    end
    if deleteplot
        uimenu(cmenu, 'label', 'delete', 'callback', @removeTab);  
    end
    if changecmap
        uimenu(cmenu, 'label', 'change colormap', 'callback', @changeColormap);  
    end
    
    %set the context menu to be used for saving each figure, must set
    %separately for all figure components, even objects within axes
    %like plot lines
    allComponents = findobj(thisTab);
    set(allComponents, 'uicontextmenu', cmenu);
        
    function saveTabbedPlot(hObj, evtData)
        %saves the figure in the current tab: creates a new figure, moves
        %the axes to it, saves, then moves the axes back to the tabs and
        %closes the figure.
        fileTypes = {'*.fig'; '*.jpg'; '*.tiff'; '*.pdf'; '*.png'}; 
        [fileName filePath filterIndex] = uiputfile(fileTypes);  
        if filterIndex ~= 0
            tabs = (get(hTabGroup, 'children'));
            selectedTab = tabs(get(hTabGroup, 'SelectedIndex'));
            cmap = get(gcf, 'colormap');
            saveFig = figure('visible', 'off');
            selectedAxes = get(selectedTab, 'children');
            set(selectedAxes, 'parent', saveFig); %move to new figure
            colormap(cmap)
            saveas(saveFig, [filePath fileName]); %save figure      
            set(selectedAxes, 'parent', selectedTab); %move back
            close(saveFig)
            %saveas function clears contextmenu, so must reset
            AddTabContextMenu(selectedTab);
        end
    end

    function removeTab(hObj, evtData)
        tabs = (get(hTabGroup, 'children'));
        selectedTab = tabs(get(hTabGroup, 'SelectedIndex'));
        delete(selectedTab);
    end


    function changeColormap(hObj, evtData)
        tabs = (get(hTabGroup, 'children'));
        selectedTab = tabs(get(hTabGroup, 'SelectedIndex'));
        selectedAxes = get(selectedTab, 'children');
        cmapList = {'jet', 'bone', 'hsv', 'hot', 'cool', 'gray'};
%         cmapIdx = find(strcmp(cmapList, inputs.colormap));
%         [cmapIdx, ok] = listdlg('liststring', cmapList, 'selectionmode', 'single', 'listsize', [150 105], 'name', 'Select colomap');
        [cmapIdx, invCmap] = ColormapDlg('liststring', cmapList,'name', 'Select Colormap');  
        if ~isempty(cmapIdx)
%             unfreezeColors
            if invCmap
                colormap(flipud(eval(cmapList{cmapIdx})));
            else
                colormap(cmapList{cmapIdx});
            end
%             for ax = selectedAxes'
%                 if strcmp(get(ax, 'tag'), 'Colorbar')
%                      cbfreeze(ax);
%                 else
%                    freezeColors(ax);              
%                 end
%             end     
        end
    end
end