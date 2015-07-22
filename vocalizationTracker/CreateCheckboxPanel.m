function pHandles = CreateCheckboxPanel(panelHandle, inputs)


    boxLength = 175;
    column1 = 10;
    column2 = 220;
    column3 = 395;
    bcolor = [0.8 0.8 0.8];
    
    %tick boxes
    pHandles.synSpecH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column1 10 boxLength 20], 'string', 'Synthetic Spectrogram','BackgroundColor', bcolor, 'value', inputs.bools.synSpec);
    pHandles.synSpecHarmH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column1 40 200 20], 'string', 'Synthetic Spectrogram Harmonics','BackgroundColor', bcolor, 'value', inputs.bools.synSpecHarm);
    pHandles.signalCompH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column1 70 boxLength 20], 'string', 'Signal Compare','BackgroundColor', bcolor, 'value', inputs.bools.signalComp);
    pHandles.ampPhaseH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column2 10 boxLength 20], 'string', 'Amplitudes and Phases','BackgroundColor', bcolor, 'value', inputs.bools.ampPhase);
    pHandles.sigSpecH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column2 40 boxLength 20], 'string', 'Signal Spectrogram','BackgroundColor', bcolor, 'value', inputs.bools.sigSpec);
    pHandles.sigSpecHarmH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column2 70 200 20], 'string', 'Signal Spectrogram Harmonics','BackgroundColor', bcolor, 'value', inputs.bools.sigSpecHarm);
    pHandles.plotResH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column3 40 boxLength 20], 'string', 'Residuals','BackgroundColor', bcolor, 'value', inputs.bools.plotRes);
%     pHandles.invertColorH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column3 70 boxLength 20], 'string', 'Invert Colormap','BackgroundColor', bcolor, 'value', inputs.bools.invertColor);

end