function [] = BuildRelease(bda);
%BuildRelease: Converts all of the .m files to .p files
%
%   [] = BuildRelease(bd)
%
%   bd   String specifying where the root directory of the BSP toolbox
%        is located
%
%   This function automatically finds all of the toolbox subdirectories
%   in bd and compiles all of the .m files in those subdirectories to 
%   .p files. Note that the .m files are not removed. This must be 
%   completed manually.
%
%   Example: Convert all of the .m files in the subdirectory 
%   'C:/Release' to .p files.
%
%      BuildRelease('C:/Release');
%
%   Version 1.00 JM
%
%   See also pcode.

%====================================================================
% Process function arguments
%====================================================================
bd = pwd;
if exist('bda') & ~isempty(bda),
    bd = bda;
    end;    
       
%====================================================================
% Main Routine
%====================================================================    
sd = pwd;                 % Starting directory
cd(bd);                   % Move to specified directory

%====================================================================
% Error Handling
%====================================================================    
fn = which(mfilename);
id = max(find(fn==char(92)));
br = fn(1:id-1);
pd = pwd;
if strcmpi(br,pd),
    error('You almost destroyed the BSP toolbox. Do not run BuildRelease in the BSP root directory!');
    end;

%====================================================================
% Main Routine
%====================================================================    
dl = dir('BuildRelease.m');
if ~isempty(dl),
    delete('BuildRelease.m'); % Do not include this function in the release
    end;
dl = dir;                 % List everything in directory
dl = dl(3:length(dl));    % Eliminate parent and current directories from list
id = find([dl.isdir]);    % Determine which are subdirectories
dl = dl(id);              % Eliminate everything that is not a subdirectory      
nd = length(dl);          % No. of subdirectories
rd = cd;                  % Save root directory
for c1 = 1:nd,            % Main loop
    cd(dl(c1).name);      % Move to subdirectory
    if strcmpi(dl(c1).name,'Documentation'), % Skip Documentation subdirectory
        delete('*.doc');  % Remove all MS Word files - only use pdf files for documentation
        cd(rd);
        continue;
        end;    
    pcode('*.m');         % Convert all .m files to .p files
    if exist('Contents.p'),
        delete('Contents.p');
        end;
    fl = dir('*.m');
    for c2 = 1:length(fl),
        fn = fl(c2).name;
        st = sprintf('rename %s Temp.m',fn);
        dos(st);
        fr = fopen('Temp.m','r');
        fw = fopen(fn,'w');
        if fr==-1 | fw==-1,
            error(sprintf('Could not open specified files ''%s'' and temp.m.',fn));
            end;
        ln = fgetl(fr); % Get the first line of the file        
        fprintf(fw,'%s\n',ln);
        ln = fgetl(fr); % Get the second line of the file
        while length(ln)>0 & ln(1)=='%', % If it's a comment line, write it
            fprintf(fw,'%s\n',ln);
            ln = fgetl(fr); % Get the next line of the file
            end;  
        fclose(fr);
        fclose(fw);
        delete('Temp.m');
        end;
        
    cd(rd);             % Move back to the root directory
    end;
cd(sd);                 % Move back to starting directory    

