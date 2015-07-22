function [varargout] = Sequence(varargin)
%Sequence: Sequence method baroreflex sensitivity (BRS).
%
%   [brs,bre,time,type,RAMPS] = Sequence(abp,abp_fs,sbpInd,ecg_fs, ...
%                            rrInd,sbpDelta,rrDelta,lag,pt,pf);
%
%   abp        Arterial blood pressure signal (mmHg).
%   abp_fs     Sampling rate of blood pressure (Hz).
%   sbpInd     Indices of systolic blood pressure values (samples).
%   ecg_fs     Sampling rate of ECG (Hz).
%   RRInd      Indices of RR values in ECG signal (samples).
%   sbpDelta   Minumum value that the systolic blood pressure must 
%              change in order for sequential values to qualify as
%              elements in a SBP ramp (mmHg). default = 1 
%              mmHg.
%   rrDelta    Minumum value that the RR interval must change in 
%              order for sequential intervals to qualify as
%              elements in a RR interval ramp (s). default = 5 ms.
%   lag        A scalar representing the number of beats after the
%              beginning of a SBP ramp to start looking for a 
%              RRI ramp (beats). Common values are [0,1,2].
%              default = 1.
%   pt         Plot type default = 0 (1 = include SBP ramps, 0 = 
%              do not include SBP ramps)
%   pf         Plot flag: 0=none (default), 1=screen.
%
%   brs    Baroreflex sensitivity index
%   bre    Baroreflex effectiveness
%   time   Time of BRS estimate, spontaneous SBP/RRI sequence 
%   type   Type = up or down sequence
%   RAMPS  SBP ramps matrix [ind length type slope time]
%          where ind is a vector of indices in the detected
%          systolic blood pressure vector not the arterial
%          blood pressure vector given
%
%   This function calculates the baroreflex sensitivity index of
%   a dual blood pressure and ECG signal using the sequence technique.
%   The sequence technique searches for spontaneous ramps in systolic
%   blood pressure and once found evaluates the RR intervals just
%   after the beginning of the ramp for a suspected baroreflex 
%   response in the form of decreaed heart rate. Spontaneous 
%   systolic blood pressure drops are evaluated in a simular 
%   fashion.
%
%   The indices of the detected systolic blood pressure peaks should
%   be passed in the vector sbpInd and in a simular manner the indices
%   of the R peak in the electrocardiogram signal should also be
%   passed in the vector rrInd.
%
%   The signals must begin at the same time and be the same length
%   for the Sequence function to return a BRS estimate that is
%   significant. The default values of the sbpDelta and rrDelta 
%   are as suggested in [Di Rienzo] for a human. However for animals
%   or infants these defalut values are not valid and should be
%   adjusted accordingly.
%
%   Example: Calculate BRS and plot results.
%
%      load dualECGABP;
%      sbpInd = PressureDetectInterbeat(abp,abp_fs);
%      rrInd  = ECGDetectRInterbeat(ecg,ecg_fs,ecg_fs);
%      brs    = Sequence(abp,abp_fs,sbpInd,ecg_fs,rrInd, ...
%                        0.1,0.002,1,1,1);
%
%   Marco Di Rienzo, Gianfranco Parati, Paolo Castiglioni, Roberto 
%   Tordi, Giuseppe Mancia, and Antonio Pedotti, "Baroreflex 
%   effectiveness index: an additional measure of baroreflex control 
%   of heart rate in daily life" Am J Physiol Regulatory Integrative 
%   Comp Physiol, vol. 280: pp. R744-R751, 2001.
%
%   Version 0.01.05.21 DT
%
%   See also ECGDetectRInterbeat and PressureDetectInterbeat.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argument Extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 5 | nargin > 10 | nargout > 5
        help('Sequence')
        return
    end
    
    % extract input arguments
    abp    = varargin{1};
    abp_fs = varargin{2};
    sbpInd = varargin{3};
    ecg_fs = varargin{4};
    rrInd  = varargin{5};


    if nargin >= 6
        sbpDelta = varargin{6};
    else
        sbpDelta = 1;
    end
    
    if nargin >= 7
        rrDelta = varargin{7};
    else
        rrDelta = 0.005;
    end
    
    if nargin >= 8
        lag = varargin{8};
    else
        lag = 1;
    end
    
    if nargin >= 9
        pt = varargin{9};
    else
        pt = 0;
    end
    
    if nargin == 10
        pf = varargin{10};
    else
        if nargout == 0
            pf = 1;
        else
            pf = 0;
        end    
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error checking
% need to expand for all input args
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sigSizeAbp = size(abp);           % sizes of input vectors
    sbpSize    = size(sbpInd);
    rrSize     = size(rrInd);

    
    if sigSizeAbp(1) == 0
        error('abp is empty.');
        return
    elseif sigSizeAbp(2) ~= 1
        error('abp is not a column vector.');
        return
    elseif isa(abp,'numeric') ~= 1   
        error('abp must be numeric');
        return
    elseif length(abp_fs)~=1
        error('abp_fs must be scalar.');
        return
    elseif isa(abp_fs,'numeric') ~= 1 
        error('abp_fs must be numeric');
        return
    elseif abp_fs <= 0
        error('abp_fs must be positive.');
        return
    elseif floor(abp_fs) ~= abp_fs
        error('abp_fs must be an integer.');
        return
    elseif length(ecg_fs)~=1
        error('ecg_fs must be scalar.');
        return
    elseif isa(ecg_fs,'numeric') ~= 1 
        error('ecg_fs must be numeric');
        return
    elseif ecg_fs <= 0
        error('ecg_fs must be positive.');
        return
    elseif floor(ecg_fs) ~= ecg_fs
        error('ecg_fs must be an integer.');
        return
    elseif sbpSize(1) == 0
        error('SBPInd is empty.');
        return
    elseif sbpSize(2) ~= 1
        error('sbpInd is not a column vector.');
        return
    elseif sbpSize(1) > sigSizeAbp(1)
        error('sbpInd is larger than abp');
        return
    elseif isa(sbpInd,'numeric') ~= 1 
        error('sbpInd must be numeric');
        return
    elseif rrSize(1) == 0
        error('rrInd is empty.');
        return
    elseif rrSize(2) ~= 1
        error('rrInd is not a column vector.');
        return
    elseif rrSize(1) > sigSizeAbp(1)
        error('rrInd is larger than abp');
        return
    elseif isa(rrInd,'numeric') ~= 1   
        error('rrInd must be numeric');
        return
    elseif length(rrDelta)~=1
        error('rrDelta must be scalar.');
        return
    elseif isa(rrDelta,'numeric') ~= 1 
        error('rrDelta must be numeric');
        return
    elseif rrDelta < 0
        error('rrDelta must not be negative.');
        return
    elseif length(sbpDelta)~=1
        error('sbpDelta must be scalar.');
        return
    elseif isa(sbpDelta,'numeric') ~= 1 
        error('sbpDelta must be numeric');
        return
    elseif sbpDelta < 0
        error('sbpDelta must not be negative.');
        return
    elseif length(lag)~=1
        error('lag must be scalar.');
        return
    elseif isa(lag,'numeric') ~= 1 
        error('lag must be numeric');
        return
     elseif floor(lag) ~= lag 
        error('lag must be an integer');
        return   
    elseif lag ~= 0 & lag ~= 1 & lag ~= 2
        warning('lag out of recommended range!');
    elseif length(pf) ~= 1  
        error('Error: pf must be scalar.');
        return
    elseif isa(pf, 'numeric') ~= 1 
        error('pf must be numeric');
        return 
    elseif pf ~= 0 & pf ~= 1
        error('pf must be a Boolean Value.');
        return
    elseif 1/ecg_fs > rrDelta
        error('rrDelta cannot be greater than time resolution.');
        return
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find SBP ramps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rampLenU = 0;                  % ramp length UP
    rampLenD = 0;                  % ramp length DOWN
    numSBP   = length(sbpInd);     % number of SBP values
    RAMPS    = [];                 % ramp [indices; lengths; codes; slopes; time] 
    UP_      = 1;                  % ramp codes 1 = UP, 2 = down
    DOWN_    = 2;
    
    for k = 2:numSBP
        
        % up ramps
        if abp(sbpInd(k)) > abp(sbpInd(k-1)) + sbpDelta
            rampLenU = rampLenU+1;
        elseif rampLenU >= 4
            regCoef  = [ ones(rampLenU,1) ...
                         sbpInd(k-rampLenU:k-1)/abp_fs] \ ...
                         abp(sbpInd(k-rampLenU:k-1));
            RAMPS = [ RAMPS [k-rampLenU;rampLenU;UP_;regCoef(2); ...
                             ((sbpInd(k-rampLenU)+sbpInd(k-1))-1)/(2*abp_fs)] ];
            rampLenU = 0;
        else
            rampLenU = 0;
        end
        
        
       % down ramps 
        if abp(sbpInd(k)) + sbpDelta < abp(sbpInd(k-1))
            rampLenD = rampLenD+1;
        elseif rampLenD >= 4
            regCoef  = [ ones(rampLenD,1) ...
                         sbpInd(k-rampLenD:k-1)/abp_fs] \ ...
                         abp(sbpInd(k-rampLenD:k-1));
            RAMPS    = [ RAMPS [k-rampLenD;rampLenD;DOWN_;regCoef(2); ...
                               ((sbpInd(k-rampLenD)+sbpInd(k-1))-1)/(2*abp_fs)] ];
            rampLenD = 0;
        else
            rampLenD = 0;
        end
    
    end
    
    % beginings of ramps in BP signal
    if length(RAMPS) == 0
        warning('Warning: No SBP ramps found');
    end
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate some statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rmSize  = size(RAMPS);          % size of ramp matrix
    num     = rmSize(2);            % number of ramps detected
    BRS     = [];                   % empty Baroreflex sensitivity
    bre     = 0;                    % Baroreflex effectiveness index
    rPtr    = 1;                    % pointer to ecg signal
    rampLen = 0;                    % length of a SBP ramp
    
    for k = 1:num
        
        startInd  = RAMPS(1,k);                   % start index in SBP signal
        startTime = (sbpInd(startInd)-1)/abp_fs;  % starting time
        rampLen   = RAMPS(2,k);                   % number of ramp values
        rampType  = RAMPS(3,k);                   % UP_ or DOWN_
        
        
        % find ECG interval
        while rPtr < length(rrInd)
            if (rrInd(rPtr)-1)/ecg_fs > startTime;
                rPtr = rPtr+lag+1;
                break;
            end
            rPtr = rPtr+1;
        end
        
        % increment to following R peak
        abpRamp = abp(sbpInd(startInd:startInd+rampLen-1));
        
        if rPtr+rampLen-1 <= length(rrInd)
          
            % find suspected RR ramp values
            rrRamp = (rrInd(rPtr:rPtr+rampLen-1) - ...
                      rrInd(rPtr-1:rPtr+rampLen-2))/ecg_fs;
                      
            % is this an RR ramp
            for ind = 2:rampLen
                if rampType == UP_
                    if rrRamp(ind-1) + rrDelta > rrRamp(ind)
                        rampLen = ind-1;
                        break;
                    end
                elseif rampType == DOWN_
                    if rrRamp(ind-1) - rrDelta < rrRamp(ind) 
                        rampLen = ind-1;
                        break;
                    end
                end    
            end 

            if rampLen >= 4
            
                regCoef  = [ ones(rampLen,1) ...
                             abpRamp(1:rampLen) ] \ ...
                             rrRamp(1:rampLen);
                % record statistics
                BRS = [ BRS ; [ regCoef(2)  ...    % BRS
                        ((sbpInd(RAMPS(1,k)) + ...  % time 
                         sbpInd(RAMPS(1,k)+rampLen-1))-1) / ...
                         (2*abp_fs) ] ...
                         rampType ];               % type

            end
            
        end   
    end
    
    if length(BRS) == 0
        warning('Warning: no RRI ramps found');
        BRS = [ 0 0 0 0 ];
    else
        bre = length(BRS(1,:))/num;               % baroreflex effectiveness
    end
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    if pf
    
        figure
        
        if pt == 1
            numplots = 4;
        else
            numplots = 3;
        end
        
        % SPB Plot
        subplot(numplots,1,1)
            plot((sbpInd-1)/abp_fs,abp(sbpInd),'.y')
            set(gca,'XTickLabel',[]);
            title('Systolic Blood Pressure')
            ylabel('SBP (mmHg)');
    
        % RR interval plot
        subplot(numplots,1,2)
            plot((rrInd-1)/ecg_fs,[rrInd(1) ; [(rrInd(2:length(rrInd))- ... 
                        rrInd(1:length(rrInd)-1))]]./ecg_fs,'.g');
            set(gca,'XTickLabel',[]);
            title('RR Interval')
            ylabel('RR Interval (s)')
            
        % BRS plot
        subplot(numplots,1,3)
        
            % indices of sequences/ramps
            upSeqInd = find( BRS(:,3) == UP_ );
            dnSeqInd = find( BRS(:,3) == DOWN_ );
            
            
            plot(BRS(upSeqInd,2),1000*BRS(upSeqInd,1),'or', ...
                 BRS(dnSeqInd,2),1000*BRS(dnSeqInd,1),'*b');
                        titStr = sprintf('Sequence method baroreflex sensitivity: BRE = %f',bre);
            title(titStr)
            ylabel('BRS (ms/mmHg)');
            legend('Up seq','Down seq');
            set(gca,'XLim',[0 length(abp)/abp_fs]);
            
            if pt == 1
                upRampInd = find( RAMPS(3,:) == UP_ );
                dnRampInd = find( RAMPS(3,:) == DOWN_ );
                subplot(numplots,1,4);
                plot(RAMPS(5,upRampInd),RAMPS(4,upRampInd),'.r', ...
                     RAMPS(5,dnRampInd),abs(RAMPS(4,dnRampInd)),'.b');
                ylabel('Slope mmHg/s');
                legend('Up ramp','Down ramp');
            end

            xlabel('Time (s)');

            FigureSet(1);
            
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargout >= 1    
        varargout{1} = BRS(:,1);    
    end
    
    if nargout >= 2
        varargout{2} = bre;
    end
    
    if nargout >= 3    
        varargout{3} = BRS(:,2);    
    end
    
    if nargout >= 4    
        varargout{4} = BRS(:,3);
    end
    
    if nargout == 5
        varargout{5} = RAMPS;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%