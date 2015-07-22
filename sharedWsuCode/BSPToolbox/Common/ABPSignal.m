function [yn,y,r] = ABPSignal(fsa,Na,fra,fca,A0a,A1a,A2a,ua,sa,sna,cfa,pfa);
%ABP: Generates a synthentic arterial blood pressure signal
%
%   [ABPLfn,ABPL,ABPg,r] = ABP(fs,N,fr,fc,A0,A1,A2,u,s,sn,cf,pf);
%   
%   fs   Signal sample rate (Hz). Default=125 Hz      
%   N    Signal duration (s). Dafault = 10 s 
%   fr   Respiratory rate, Hz. Default = 0.4 Hz
%   fc   Heart rate (Hz). Default = 1.9 Hz
%   A0   Respiratory amplitude (0-1). Default = 0.1 
%   A1   Incident wave amplitude (0-1). Default = 1 
%   A2   Reflected wave amplitude (0-1). Dafault = 0.5 
%   u    Mean arterial pressure (mmHg). Default = 80 mmHg
%   s    Arterial pressure standard deviation (mmHg). Default = 15 mmHg
%   sn   Noise standard deviation (mmHg). Default = s/2 mmHg
%   cf   LTI filter cutoff frequency (Hz). Default = 10 Hz
%   pf   Plot flag: 0=none (default), 1=screen
%
%   yn   ABP signal with noise
%   y    ABP signal without noise
%   r    Respiration
%
%   Creates a synthetic arterial blood pressure signal    
%
%   Example: Create 20 seconds of an arterial blood pressure signal sampled 
%   at 100 Hz, with the cardiac component at 2 Hz and the respiratory component 
%   at 0.5 Hz. 
% 
%      ABPSignal(100,20,0.5,2);
%
%   Version 1.00 MA
%
%   See also TBC

%=====================================================================
% Process function arguments
%=====================================================================
if nargin<1 | nargin>12,
    help ABPSignal;
    return;
    end;
    
fs    = 125;                          % Sampling Frequency, Hz
T     = 1/fs;                         % Sampling Perid, s  
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    T  = 1/fs;
end;  
     
N     = 10*fs;                        % Signal Duration, s
n     = 0:N-1;                        % Index  
if exist('Na') & ~isempty(Na),
    N     = Na*fsa;                   
    n     = 0:N-1;                    
end;

fr    = 0.5;                          % Repiratory Rate, Hz
if exist('fra') & ~isempty(fra),
    fr = fra; 
end;

fc    = 1.9;                          % Heart Rate (cardiac frequency), Hz
if exist('fca') & ~isempty(fca),
    fc = fca; 
end;

A0    = 0.1;                          % Respiration Amplitude
if exist('A0a') & ~isempty(A0a),
    A0 = A0a; 
end;

A1    = 1;                            % Incident Wave Amplitude
if exist('A1a') & ~isempty(A1a),
    A1 = A1a; 
end;

A2    = 0.7;                          % Reflected Wave Amplitude 
if exist('A2a') & ~isempty(A2a),
    A2 = A2a; 
end;

u     = 80;                           % Mean Arterial Pressure, mmHg
if exist('ua') & ~isempty(ua),
    u = ua; 
end;

s     = 15;                           % Standard Deviation, mmHg 
if exist('sa') & ~isempty(sa),
    s = sa; 
end;

sn    = s/10;                          % Standard Deviation (Noise)   
if exist('sna') & ~isempty(sna),
    sn = sna; 
end;

cf    = 7;                        % LTI Filter cutoff frequency 
if exist('cfa') & ~isempty(cfa),
    cf = cfa; 
end;

pf = 0;                                % Default - no plotting
if nargout==0,                         % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%==========================================================================
% User Specified Parameters (script)
%==========================================================================
%fs    = 125;                          % Sampling Frequency, Hz
%T     = 1/fs;                         % Sampling Perid, s       
%N     = 10*fs;                        % Signal Duration, s
%n     = 0:N-1;                        % Index  

%fr    = 0.4;                          % Repiratory Rate, Hz
%fc    = 1.9;                          % Heart Rate (cardiac frequency), Hz
%A0    = 0.1;                          % Respiration Amplitude
%A1    = 1;                            % Incident Wave Amplitude
%A2    = 0.4;                          % Reflected Wave Amplitude 
%u     = 80;                           % Mean Arterial Pressure, mmHg
%s     = 15;                           % Standard Deviation, mmHg 
%sn    = s/2;                          % Standard Deviation (Noise)   
%fc    = 10;                           % LTI Filter cutoff frequency 
         

%==========================================================================
% Waveform Generation
%==========================================================================
r     = A0*cos(2*pi*fr.*n*T);         % Respiration
iw    = A1*cos(2*pi*fc.*n*T);         % Incident Wave  
rw    = A2*cos(2*pi*2*fc.*n*T+pi/2);    % Reflected Wave
ABP   = [r+10*abs(min(r))].*[iw+rw];  % ABP signal (locally)
ABPn  = (ABP-mean(ABP))./std(ABP);    % ABP signal normalized
ABPL  = ABPn*s+u;                     % ABP(u,s;t) signal (locally)

noise  = sn*randn(1,N);               % Gaussian Noise
ABPLg  = ABPL+noise;                  % Additive Noise  
ABPLfn = lowpass(ABPLg,fs,cf);        % ABP signal + noise

yn = ABPLfn+r;
y  = ABPL+r;

%==========================================================================
% Plotting
%==========================================================================
plot(yn)
spectrogram(yn,20,blackman(round(10*fs)));

pfp = 0;
if pfp==1
    figure; figureset(1)
    h = plot(n/fs,ABPLfn);
    xlabel('Time,s');
    ylabel('ABP, mmHg');
    box off; axisset;
    
    figure; figureset(2);
    h = plot(n/fs,ABPLfn,'b', n/fs, ABPL,'r');
    xlabel('Time,s');
    ylabel('ABP, mmHg');
    box off; axisset;
    
        
    figure; figureset(4);
    h = plot(n/fs,iw,'b', n/fs, rw,'r', n/fs, r);
    xlabel('Time,s');
    ylabel('ABP Components, normalized');
    box off; axisset(8); ylim([-1.2 1.2]);
    
    figure; figureset(1,'wide');
    subplot(3,1,1);
    h = plot(n/fs,iw,'b', n/fs, rw,'r', n/fs, r);
    xlabel('Time,s');
    ylabel('ABP Components, normalized');
    box off; axisset(8); ylim([-1.2 1.2]);
    
    subplot(3,1,2)
    h = plot(n/fs, ABPL);
    ylabel('ABP_L Components, mmHg');
    box off; axisset(8);
    
    subplot(3,1,3)
    h = plot(n/fs,ABPLfn);
    xlabel('Time,s');
    ylabel('ABP_{LN}, mmHg');
    box off; axisset(8);
end