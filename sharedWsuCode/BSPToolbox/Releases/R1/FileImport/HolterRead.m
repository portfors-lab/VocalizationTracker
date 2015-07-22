function [y,fs] = QRSRead(fn,ld,npa,pf);
%HolterRead: Reads Holter ECG data.
%
%   [y,fs] = QRSRead(fn,ld,np,pf);
%
%   fn   File name (string).
%   ld   Lead (integer, 1-3). Default = 1.
%   np   Number of points to read. Default = inf.
%        If a vector, specifies first and last point to read.
%   pf   Plot flag. 
%
%   y    Requested signal segment of the electrocardiogram (scaled).
%   fs   Sample rate (Hz).
%
%   Reads the Holter waveform data for the filename specified. 
%   Returns a matrix with three columns where each column contains 
%   the QRS waveform from each of the three leads. The sample rate is 
%   fixed at 128 Hz.
% 
%   This was developed based on data acquired from a Holter monitor
%   and some careful reverse-engineering of the file format. It may 
%   not work with all data acquired from Holter monitors.
%
%   If the size is specified, this only reads in the first Size 
%   samples. Otherwise, it reads the entire file.
% 
%   If PlotFlag is included, will generate plot of each lead.
