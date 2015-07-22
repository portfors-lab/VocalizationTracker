try
    savestr = [workspaceName '_' num2str(now) '.mat'];
catch
    display('workspaceName variable not specified for');
    display('saveWorkspace script. Defaulting to ''workspace''');
    workspaceName = 'workspace';
    savestr = [workspaceName '_' num2str(now) '.mat'];
end
if exist(savestr,'file')
    savestr = [savestr(1:end-4) '_b.mat'];
end
save(savestr);
display(['Workspace saved as ' savestr]);