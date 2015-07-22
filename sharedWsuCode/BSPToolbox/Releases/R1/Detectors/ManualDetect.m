function [mp,ibi] = ManualDetect(x);
%ManualDetect: Manual Peak Detector 
%   
%   [mp,ibi] = ManualDetect(x)
%
%   x     Input signal            
%
%   mp    Manual Detected Peaks, samples
%   ibi   Interbeat interval, samples
%
%   ManualDetect can be used to perform manual annotation of 
%   wave-component on signals such as ICP, ABP, ECG. The
%   function calculates the interbeat intervals, the first, 
%   and the second derivative, and plots them to aid
%   the manual annotation in regions with high artifact. 
%
%   Example: Perform manual detection on an the ICP signal. 
%
%      load ICP; 
%      [mp,ibi] = ManualDetect(icp);
%
%   Version  1.00 MA
%
%   See also PressureDetect, ECGDetectQRS and ECGDetectRInterbeat. 
