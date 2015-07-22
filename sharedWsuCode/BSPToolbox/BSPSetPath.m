function [] = BSPSetPath(bd);
%BSPSetPath: Adds the BSP toolbox to MATLAB's path
%
%   [] = BSPSetPath(bd)
%
%   bd   String specifying where the root directory of the BSP toolbox
%        is located
%
%   This should be added to the startup.m file unique to each user that
%   MATLAB calls when it first starts up.
%
%   Example: Add the BSP toolbox to MATLAB's path. The user has placed
%   the BSP toolbox in the directory 'C:/BSPToolbox'.
%      BSPSetPath('C:/BSPToolbox');
%
%   Version 1.00 JM
%
%   See also pwd, addpath, and path.

addpath([bd '/Common'],'-end'); 
addpath([bd '/Correlation'],'-end');
addpath([bd '/Data'],'-end');
addpath([bd '/Detectors'],'-end');
addpath([bd '/Documentation'],'-end'); 
addpath([bd '/FileImport'],'-end'); 
addpath([bd '/Metrics'],'-end'); 
addpath([bd '/Models'],'-end'); 
addpath([bd '/NonlinearFilters'],'-end');
addpath([bd '/SpectralAnalysis'],'-end');
addpath([bd '/Visualization'],'-end');
addpath([bd '/FileImport/WindowsBinaries'],'-end');


