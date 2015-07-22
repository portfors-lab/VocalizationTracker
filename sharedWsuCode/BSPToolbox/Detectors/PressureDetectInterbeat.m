function [B,A,C,D] = PressureDetectInterbeat(x,fsa,wla,pfa);
%PressureDetectInterbeat: Pressure detector algorithm (ICP & ABP)
%
%   [B, A, C, D] = PressureDetectInterbeat(x,fs,wl,pf);
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=125 Hz      
%   wl   Window length (s), default=5 s 
%   pf   Plot flag: 0=none (default), 1=screen
%
%   B    Percusion peak (P1) index, samples
%   A    Minima preceding the percusion peak (P1) index, samples
%   C    Dichrotic notch (DN) index, samples
%   D    Dichrotic peak (P2) index, samples
%
%   Beat detection algorithm for pressure signals based on ranks and 
%   intebeat distance. The algorithm detects the minima prior to the 
%   percursor peak (A), the percursor peak (B), the dicrotic notch (C),
%   and the dicrotic peak (D). The algorithm (1) detects all maxima and
%   minima in the raw data, (2) detrends and smooths the data using a 
%   bandpass filter, (3) estimates the heart rate of the input signal as
%   a function of time (window), (4) detects and classifies all maxima 
%   and minima in each window of filtered signal based on ranks, and 
%   (5) corrects the missed beats and false positives based on the 
%   interbeat distance. 
%
%   Example: Find the minima (A) prior to the percursor peak, the percursor
%   peak (B), the dicrotic notch (C), and the dicrotic peak (D) in the ICP 
%   signal
%
%      load ICPCOMPOSITE; 
%      PressureDetectInterbeat(icpcomposite, fs);
%
%   Version 1.00 MA
%
%   See also PressureDetect, PressureDetectRank, and ECGDetectQRS.


%===========================================================================
% Process function arguments
%===========================================================================
if nargin<1 | nargin>4,
    help PressureDetectInterbeat;
    return;
end;

fs = 125;                                        % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
end;

wl = 5*fs;                                      % Default sampling rate, Hz
if exist('wla') & ~isempty(wla),
    wl = wla*fs;
end;
        
pf = 0;                                          % Default - no plotting
if nargout==0,                                           % Plot if no output arguments
    pf = 1;
end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
end;

x  = x(:);
LD = length(x);

%===========================================================================
% Detect maxima and minima on raw data
%===========================================================================
MN     = DetectMinima(x);
PK     = DetectMaxima(x);


%===========================================================================
% Detrend and smooth data
%===========================================================================
ysm    = lowpass(x, fs, 30, 1, 2, 0);   % Filter out all frequencies > 30
yntr   = highpass(ysm, fs, 0.4, 2, 0);  % Filter out all frequencies < 0.4


%===========================================================================
% Estimate HR and T_hat
%===========================================================================
[hr] = HeartRateSpectral(x,fs,10,0);
mhr    = median(hr);
T_est  = (1./hr)*fs;
T_hat  = (1/mhr)*fs; 


%===========================================================================
% Detection & Classification based on ranks
%===========================================================================
A = [];
B = [];
C = [];
D = [];

for k=1:wl:LD
    if ((k+wl-1 < LD-wl))
        yntrs= yntr(k:k+wl-1);
        
        thra = prctile(yntrs, 10);
        As   = DetectMinima(yntrs, thra);    
        As   = As(:);
        A    = [A; As+k-1];
        
        thrb = prctile(yntrs, 80);
        Bs   = DetectMaxima(yntrs, thrb);
        Bs   = Bs(:);   
        B    = [B; Bs+k-1];
   
        Cs   = DetectMinima(yntrs);
        Cs   = Cs(:);
        C    = [C; Cs+k-1];
        
        Ds   = DetectMaxima(yntrs);
        Ds   = Ds(:);
        D    = [D; Ds+k-1];

    else
        yntrs = yntr(k:LD);
        
        thra = prctile(yntrs, 10);
        As   = DetectMinima(yntrs, thra);    
        As   = As(:);
        A    = [A; As+k-1];

        thrb = prctile(yntrs, 80);
        Bs   = DetectMaxima(yntrs, thrb);
        Bs   = Bs(:);   
        B    = [B; Bs+k-1];
               

        Cs   = DetectMinima(yntrs);
        Cs   = Cs(:);
        C    = [C; Cs+k-1];
       
        Ds   = DetectMaxima(yntrs);
        Ds   = Ds(:);
        D    = [D; Ds+k-1];
    end;
end;

Aa = unique(A);
Ba = unique(B);
Ca = unique(C);
Da = unique(D);
Aa = sort(Aa);
Ba = sort(Ba);
Ca = sort(Ca);
Da = sort(Da);


%===========================================================================
% Correction based on the interbeat interval 
%===========================================================================

%---------------------------------------------------------------------------
% Correction for A
%---------------------------------------------------------------------------
A2Aa = diff(Aa);
A2A  = A2Aa; 
Am   = mean(A2Aa);
T    = median([T_est; Am]);
LL   = (1/2)*T_hat;

%..........................................................................
% Correct false positives
%..........................................................................
I    = find(A2Aa < LL);
I    = I+1;
LI   = length(A);
if ~isempty(I),   
    for k=1:LI,
        I    = find(A~=A(I(1)));
        A    = A(I);
        A    = unique(A);
        A    = sort(A);
        A2A  = diff(A);
        T    = mean(A2A);
        T    = median([T; T_est]);
        LL   = (1/2)*T;
        I    = find(A2A < LL);
        I    = I+1;
        if isempty(I)
            break
        end
       
    end;
end;


%............................................................................
% Correct misdetections
%............................................................................
UL   = (2)*T;
I    = find(A2A > UL);
if ~isempty(I),
    for k=1:LI,  
        ni = A(I(1))+T;
        [mag, I] = min(abs(ni-MN));
        A   = [A; MN(I)];
        A   = unique(A);
        A   = sort(A);
        A2A = diff(A);
        T   = mean(A2A);
        T   = median([T; T_est]);
        UL  = 1.5*T;
        I = find(A2A > UL);
        if isempty(I),
            break
        end;
    end;
end;



%---------------------------------------------------------------------------
% Correction for B
%---------------------------------------------------------------------------
B2Ba = diff(Ba);
B2B  = B2Ba;
T    = T;
LL   = (1/2)*T;

%..........................................................................
% Correct false positives
%..........................................................................
I    = find(B2Ba < LL);
I    = I+1;
LI   = length(B);
if ~isempty(I),   
    for k=1:LI,
        I    = find(B~=B(I(1)));
        B    = B(I);
        B    = unique(B);
        B    = sort(B);
        B2B  = diff(B);
        T    = mean(B2B);
        T    = median([T; T_est]);
        LL   = (1/2)*T;
        I    = find(B2B < LL);
        I    = I+1;
        if isempty(I)
            break
        end
       
    end;
end;

    
%............................................................................
% Correct misdetections
%............................................................................
UL   = (2)*T;
I    = find(B2B > UL);
if ~isempty(I),
    for k=1:LI,  
        ni = B(I(1))+T;
        [mag, I] = min(abs(ni-PK));
        B   = [B; PK(I)];
        B   = unique(B);
        B   = sort(B);
        B2B = diff(B);
        T   = mean(B2B);
        T   = median([T; T_est]);
        UL  = 1.5*T;
        I = find(B2B > UL);
        if isempty(I),
            break
        end;
    end;
end;


                                         
%---------------------------------------------------------------------------
% Correction for C
%---------------------------------------------------------------------------  
for i=1:length(A);
    C   = C(C~=A(i));
end;

Ca   = unique(C);
C2Ca = diff(Ca);
C2C  = C2Ca;

T    = T;
LL   = (1/2)*T;

%..........................................................................
% Correct false positives
%..........................................................................
I    = find(C2Ca < LL);
I    = I+1;
LI   = length(C);
if ~isempty(I),   
    for k=1:LI,
        I    = find(C~=C(I(1)));
        C    = C(I);
        C    = unique(C);
        C    = sort(C);
        C2C  = diff(C);
        LL   = (1/2)*T;
        I    = find(C2C < LL);
        I    = I+1;
        if isempty(I)
            break
        end
       
    end;
end;


UL   = (2)*T;

%............................................................................
% Correct misdetections
%............................................................................
I    = find(C2C > UL);
if ~isempty(I),
    for k=1:LI,  
        ni = C(I(1))+T;
        [mag, I] = min(abs(ni-MN));
        C   = [C; MN(I)];
        C   = unique(C);
        C   = sort(C);
        C2C = diff(C);
        UL  = 1.5*T;
        I = find(C2C > UL);
        if isempty(I),
            break
        end;
    end;
end;


for i=1:length(A);
    C   = C(C~=A(i));
end;


%---------------------------------------------------------------------------
% Correction for D
%--------------------------------------------------------------------------- 
for i=1:length(B);
    D   = D(D~=B(i));
end;

Da   = D;
D2Da = diff(Da);
D2D  = D2Da;

T    = T;
LL   = (1/4)*T;

%..........................................................................
% Correct false positives
%..........................................................................
I    = find(D2Da < LL);
I    = I+1;
LI   = length(D);
if ~isempty(I),   
    for k=1:LI,
        I    = find(D~=D(I(1)));
        D    = D(I);
        D    = sort(D);
        D2D  = diff(D);
        LL   = (1/2)*T;
        I    = find(D2D < LL);
        I    = I+1;
        if isempty(I)
            break
        end
       
    end;
end;


UL   = (1.5)*T;

%............................................................................
% Correct misdetections
%............................................................................
I    = find(D2D > UL);
if ~isempty(I),
    for k=1:LI,  
        ni = D(I(1))+T;
        [mag, I] = min(abs(ni-PK));
        D  = [D; PK(I)];
        D  = unique(D);
        D  = sort(D);
        UL  = 1.5*T;
        I = find(D2D > UL);
        if isempty(I),
            break
        end;
    end;
end;



for i=1:length(B);
    D   = D(D~=B(i));
end;



%====================================================================
% Interbeat Intervals 
%====================================================================
A2A = diff(A);
B2B = diff(B);
C2C = diff(C);
D2D = diff(D);

%====================================================================
% Plotting
%====================================================================

if pf == 1, 
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