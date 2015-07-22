function UpdateCVS(localPath_a)
%UpdateCVS: Update a working copy from the CVS repository
%
%   [] = UpdateCVS(localPath, cvsServer)
%
%   localPath         The directory the working copy resides in.
%
%   Updates the local working copy of a CVS module from the CVS 
%   repository at spike.ds.cecs.pdx.edu\CVS. This function relies
%   on the CVS executable being in the windows path.
%
%   Example: Update a local copy of the BSP Toolbox
%
%      UpdateCVS('F:\BSPToolbox');
%

localPath = localPath_a;

oldDir = pwd;
chdir(localPath);
!cvs -d :sspi:spike/CVS update
chdir(oldDir);


