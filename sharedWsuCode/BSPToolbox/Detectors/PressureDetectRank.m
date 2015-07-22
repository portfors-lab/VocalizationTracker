function [B,A,C,D] = PressureDetectRank(x,fsa,dta,pfa);
%PressureDetectRank: Pressure detector based on ranks
%
%   [B,A,C,D] = PressureDetectRank(x,fs,dt,pf);
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=125 Hz      
%   dt   Type: 0=perform detection on raw data, 1=detrend first (default) 
%   pf   Plot flag: 0=none (default), 1=screen,
%
%   B    Percusion peak (P1) index, samples
%   A    Minima preceding the percusion peak (P1) index, samples 
%   C    Dicrotic notch (DN) index, samples
%   D    Dicrotic peak (P2) index, samples
%
%   Beat detection algorithm for pressure signals based on ranks. The 
%   algorithm detects the minima (A) prior to the percursor peak, the 
%   percursor peak (B), the dicrotic notch (C), and the dicrotic peak (D).
%
%   Example: Find the minima (A) prior to the percursor peak, the percursor
%   peak (B), the dicrotic notch (C), and the dicrotic peak (D) in the ICP 
%   signal.
%
%      load ICP; 
%      [B,A,C,D] = PressureDetectRank(icp,fs,1,1);
%
%   Version 1.00 MA
%
%   See also PressureDetect, PressureDetectInterbeat, and ECGDetectQRS.


%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1 | nargin>4,
    help PressureDetectRank;
    return;
    end;

fs = 125;                                 % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;
    
dt = 1;  % Type, perform detection on the raw data
if exist('dta') & ~isempty(dta),
    dt = dta;
    end;
    
pf = 0; % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

x  = x(:);
LD = length(x);


%=========================================================================
% Detrend Data Type ~-0
%=========================================================================
if dt ~= 0,
    ysm = lowpass(x, fs, 30, 1, 2, 0);
    yntr = highpass(ysm, fs, 0.5, 2, 0);


%=========================================================================
% Detector
%=========================================================================
A = [];
B = [];
C = [];
D = [];


wl= 5*fs;

for k=1:wl:LD
    if ((k+wl-1 < LD-wl))
        yntrs= yntr(k:k+wl-1);
        
        thra = prctile(yntrs, 15);
        As   = DetectMinima(yntrs, thra);    
        As   = As(:);
        A    = [A; As+k-1];

        thrb = prctile(yntrs, 80);
        Bs   = DetectMaxima(yntrs, thrb);
        Bs   = Bs(:);   
        B    = [B; Bs+k-1];
               
        thrc = prctile(yntrs, 15);
        Cs   = DetectMinima(yntrs);
        Cs   = Cs(yntrs(Cs) > thrc);
        Cs   = Cs(:);
        C    = [C; Cs+k-1];
       
        thrd = prctile(yntrs, 85);
        Ds   = DetectMaxima(yntrs);
        Ds   = Ds(yntrs(Ds) < thrd);
        Ds   = Ds(:);
        D    = [D; Ds+k-1];

    else
        yntrs = yntr(k:LD);
        
        thra = prctile(yntrs, 15);
        As   = DetectMinima(yntrs, thra);    
        As   = As(:);
        A    = [A; As+k-1];

        thrb = prctile(yntrs, 80);
        Bs   = DetectMaxima(yntrs, thrb);
        Bs   = Bs(:);   
        B    = [B; Bs+k-1];
               
        thrc = prctile(yntrs, 15);
        Cs   = DetectMinima(yntrs);
        Cs   = Cs(yntrs(Cs) > thrc);
        Cs   = Cs(:);
        C    = [C; Cs+k-1];
       
        thrd = prctile(yntrs, 85);
        Ds   = DetectMaxima(yntrs);
        Ds   = Ds(yntrs(Ds) < thrd);
        Ds   = Ds(:);
        D    = [D; Ds+k-1];
    end;
end;
end;


%=========================================================================
% Type =0 =
%=========================================================================
if dt == 0,
    yntr = x;

    thra = prctile(yntr, 15);
    A    = DetectMinima(yntr, thra);    


    thrb = prctile(yntr, 80);
    B    = DetectMaxima(yntr, thrb);

               
    Cs   = DetectMinima(yntr);
    C    = Cs(yntr(Cs) > thra);

       
    Ds   = DetectMaxima(yntr);  
    D    = Ds(yntr(Ds) < thrb);
end;

A = unique(A);
B = unique(B);
C = unique(C);
D = unique(D);
A = sort(A);
B = sort(B);
C = sort(C);
D = sort(D);

% Interbeat difference
A2A  = diff(A);
B2B  = diff(B);
C2C  = diff(C);
D2D  = diff(D);

% Intebeat Indexes
A2Ai = (A(1:length(A)-1) + A(2:length(A)))/2;
B2Bi = (B(1:length(B)-1) + B(2:length(B)))/2;
C2Ci = (C(1:length(C)-1) + C(2:length(C)))/2;
D2Di = (D(1:length(D)-1) + D(2:length(D)))/2;

%=========================================================================
% Plotting
%=========================================================================
if pf == 1    
    tx   = 1:length(x);
    
    A2Ai = (A(1:length(A)-1) + A(2:length(A)))/2;
    B2Bi = (B(1:length(B)-1) + B(2:length(B)))/2;
    C2Ci = (C(1:length(C)-1) + C(2:length(C)))/2;
    D2Di = (D(1:length(D)-1) + D(2:length(D)))/2;
    
    figure (1)
    figureset(1)
    h=plot(tx, x, A, x(A), 'k.', B, x(B), 'r.', C, x(C), 'c.', D, x(D), 'g.');
    set(h, 'MarkerSize', 15);
    axis('tight');
    title('Detection');
    xlabel('Sample');
    ylabel('Magnitude');
    
    figure(2)
    figureset(2)
    subplot(5,1,1)
    h=plot(tx, x, A, x(A), 'k.', B, x(B), 'r.', C, x(C), 'c.', D, x(D), 'g.');
    set(h, 'MarkerSize', 15);
    axis('tight');
    title('Detection and interbeat intervals-A, B, C, D');
    
    subplot(5,1,2)
    plot(A2Ai, A2A);

    subplot(5,1,3)
    plot(B2Bi, B2B);

    subplot(5,1,4)
    plot(C2Ci, C2C);

    subplot(5,1,5)
    plot(D2Di, D2D);

    figure(3)
    figureset(3)
    h=plot(1:length(x), x, 'b', A, x(A), 'r.', A2Ai, (A2A-mean(A2A))./std(A2A));
    set(h, 'Markersize', 15);
    axis('tight');
    title('Detection: A-Component');
    xlabel('Sample');
    ylabel('Magnitude');
    
    figure(4)
    figureset(4)
    h=plot(1:length(x), x, 'b', B, x(B), 'r.', B2Bi, (B2B-mean(B2B))./std(B2B));
    set(h, 'Markersize', 15);
    axis('tight');
    title('Detection: B-Component');
    xlabel('Sample');
    ylabel('Magnitude');
    
    figure(1)
end;

%=====================================================================
% Take care of outputs
%=====================================================================
if nargout==0,
clear('A','B', 'C', 'D');
end;   
