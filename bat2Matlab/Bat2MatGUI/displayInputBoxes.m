function handles = displayInputBoxes(tree, mpanel)

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
     
        %now display input boxes for selected nodes
        %nodes = tree.getSelectedNodes;
       nodes = tree.getSelectionModel().getSelectionPaths()
       
        layout = GridLayout(0,3);
        vgap = 0;
        layout.setVgap(vgap);
        layout.setHgap(10);
        
        jPanel = JPanel(layout);

        for x = 1:length(nodes)
            %get the path string
            path = nodes(x).getPath.cell;
            subPathStrs = cellfun(@(p) [p.getName.char, filesep], path, 'un',0);
            pathStr = strrep([subPathStrs{:}], [filesep,filesep], filesep);
        
            handles(x).pathString = pathStr;
            %get the name of the folder
            rev = fliplr(pathStr);
            foldername = strtok(rev, '/');
            foldername = fliplr(foldername);
            
            handles(x).test = JTextField();
            handles(x).trace = JTextField();
            handles(x).threshold = JTextField();
            
            jPanel.add(JLabel([foldername ' test # :']));
            jPanel.add(JLabel([foldername ' trace # :']));
            jPanel.add(JLabel([foldername ' threshold :']));
            
            jPanel.add(handles(x).test);
            jPanel.add(handles(x).trace);
            jPanel.add(handles(x).threshold);
            
        end
        
        
        %determines the height of the panel due to how many folders are
        %selected, the number 47 was found by trial and error, it is the
        %necessary minimum height for the textfield + label        
        
        if exist('handles', 'var')
            set(mpanel, 'units', 'pixels');
            panelPos = get(mpanel, 'position');
            n = length(handles);
            heightDemand = n*(47) + 2*n*vgap;
            panelHeight = min(heightDemand, panelPos(4));
            scrollPanePos = [0 0 panelPos(3) panelHeight];
            scrollPane = JScrollPane(jPanel);   
            [jhandle mhandle] = javacomponent(scrollPane, scrollPanePos, mpanel);
        else
            handles = [];
        end

end