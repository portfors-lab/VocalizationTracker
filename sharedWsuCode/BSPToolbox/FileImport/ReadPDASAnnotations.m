function [C] = ReadAnnoations(fn);
%clear all;
%close all;
%fn = 'R:\PDAS\TBI\Patient007\Session001\Annotations.txt';

fid = fopen(fn,'r');
if fid==-1,
    error('Could not open file.');
    end;
    
sl = fgetl(fid);
dl = find(sl==',');
dl = [0,dl,length(sl)+1];

nc = length(dl)-1;
C = cell(100,nc);

for c1=1:length(dl)-1,
    C(1,c1) = cellstr(sl(dl(c1)+1:dl(c1+1)-1));
    end;

cc = 1;    
sl = fgetl(fid);    
while sl~=-1,
    dl = find(sl==',');
    dl = [0,dl,length(sl)+1];
    cc = cc + 1;
    for c1=1:length(dl)-1,
        C(cc,c1) = cellstr(sl(dl(c1)+1:dl(c1+1)-1));
        end; 
    for c1=length(dl):nc,
        C(cc,c1) = cellstr('Err');
        end;
    sl = fgetl(fid);    
    end;
C = C(1:cc,:);
    