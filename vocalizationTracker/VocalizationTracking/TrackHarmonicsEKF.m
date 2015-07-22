function Estimates = TrackHarmonicsEKF(y,sampleRate,ModelParameters,maxHarm, varargin)
%TrackHarmonicsEKF: Tracks sinusoidal harmonic amplitudes and frequency
%
%   StateEstimated =  TrackHarmonicsEKF(y,sampleRate,ModelParameters,
%                     'ParameterName',value,...);
%
% Input Parameters:
%   y                 Observed signal. Scalar.
%   sampleRate        Sample rate (Hz).
%   ModelParameters   Structure containing all of the parameters that 
%                     describe the statistical model.
%
% Optional Parameters:
%   'smoother'        Specifies whether to apply the smoothing recursions
%                     Default: True
% 
% Output
%   Estimates         Structure containing all of the estimated parameters
%                     and state by the EKF equations.
%
%   Uses the extended Kalman filter (EKF) and or the extended Kalman smoother
%   (EKS) recursions to estimate the
%   state of a sinusoidal harmonic model.
%
%   Example: To be written.
%
%   Reference: To be completed.
%
%   See also NonparametricSpectrogram.

%   Tags: Filter, Estimator

%==============================================================================
% Abbreviations
%==============================================================================
% x  State vector (time varying).
% y  Measurement vector (time varying).
% Q  State covariance matrix.
% R  Measurement covariance matrix. Scalar, in this case.
% H  Jacobian of the measurement function.
% F  Jacobian of the state transition/update function.
% Kf Kalman filter gain (filtered estimate).
% Pf State error covariance for the filtered state estimates.
% Pp State error covariance for the predicted state estimates.
% Re Innovation (prediction error) covariance matrix. Scalar, in this application.
% xp Predicted estimate of the state.
% xf Filtered  estimate of the state.

%==============================================================================
% Error Checking
%==============================================================================
if nargin<3
    help TrackHarmonics;
    return;
end

%==============================================================================
% Process Function Arguments
%==============================================================================
smoother = true; %Boolean indicating whether to apply the smoothing recursions

nMandatoryArguments = 3;
if nargin>nMandatoryArguments
    if ~isstruct(varargin{1})
        if rem(length(varargin),2)~=0, error('Optional input arguments must be in name-value pairs.'); end;
        Parameters = struct;
        for c1=1:2:length(varargin)-1
            if ~ischar(varargin{c1}), error(['Error parsing arguments: Expected property name string at argument ' num2str(c1+1)]); end        
            Parameters.(varargin{c1}) = varargin{c1+1};
        end
    else
        Parameters = varargin{1};
    end
    
    parameterNames = fieldnames(Parameters);
    for c1 = 1:length(parameterNames)
        parameterName  = parameterNames{c1};
        parameterValue = Parameters.(parameterName);
        switch lower(parameterName)
            case lower('smoother'),  smoother = parameterValue;            
            otherwise,               error(['Unrecognized property: ''' varargin{c1} '''']);
        end
    end
end

%==============================================================================
% Preprocessing
%==============================================================================
nSamples        = length(y);
%nHarmonics      = ModelParameters.nHarmonics;
nHarmonics      = maxHarm;
sampleInterval  = 1/sampleRate;
nStateVariables = 2+2*nHarmonics;
frequencyMean   = ModelParameters.frequencyMean;
frequencyRange  = ModelParameters.frequencyRange;
frequencyAlpha  = ModelParameters.frequencyCoefficient;
amplitudeAlpha  = ModelParameters.amplitudeCoefficient;

%==============================================================================
% Declare Structure to Store All of the Estimates In
%==============================================================================
StateErrorVariance = struct(...
    'phase'        ,nan(nSamples,1),...
    'frequency'    ,nan(nSamples,1),...
    'amplitudesCos',nan(nSamples,nHarmonics),...
    'amplitudesSin',nan(nSamples,nHarmonics)...
    );

State = struct(...
    'y'            ,nan(nSamples,1),...
    'error'        ,nan(nSamples,1),...
    'phase'        ,nan(nSamples,1),...
    'frequency'    ,nan(nSamples,1),...
    'amplitudesCos',nan(nSamples,nHarmonics),...
    'amplitudesSin',nan(nSamples,nHarmonics)...
    );

Estimates = struct(...
    'label','Extended Kalman Filter',...
    'Filtered' ,State,...
    'Predicted',State ...
    );

if smoother
    Estimates.Smoothed = State;
    Estimates.Smoothed.jacobianFunction = struct( 'state',cell(nSamples,1),...                               
                                                  'measurement',cell(nSamples,1));
    Estimates.Smoothed.ErrorCovariance =  struct( 'state',cell(nSamples,1),...                               
                                                  'measurement',cell(nSamples,1));   
end

Estimates.ModelParameters = ModelParameters;
Estimates.ModelParameters.sampleRate = sampleRate;

%==============================================================================
% Initialize Predicted State at Time 0
%==============================================================================
xp = [ModelParameters.StateInitial.phase;      ...
      ModelParameters.StateInitial.frequency;  ...
      ModelParameters.StateInitial.amplitudes; ...
      ModelParameters.StateInitial.amplitudes  ...
      ];

F      = eye(nStateVariables);
F(1,2) = 2*pi*sampleInterval;
F(2,2) = frequencyAlpha;

Q      = diag([ModelParameters.StateVariance.phase;...
               ModelParameters.StateVariance.frequency;...
               ModelParameters.StateVariance.amplitudes;...
               ModelParameters.StateVariance.amplitudes...
               ])*sampleInterval;

R      = ModelParameters.measurementVariance;

%Pp     = 0.01*eye(nStateVariables);
Pp     = diag([ModelParameters.StateVarianceInitial.phase;...
               ModelParameters.StateVarianceInitial.frequency;...
               ModelParameters.StateVarianceInitial.amplitudes;...
               ModelParameters.StateVarianceInitial.amplitudes...
               ])*sampleInterval;


%==============================================================================
% Main Loop
%==============================================================================
for n=1:nSamples
    
    %--------------------------------------------------------------------------
    % EKF Recursions
    %--------------------------------------------------------------------------    
    phasePredicted         = xp(1);
    frequencyPredicted     = frequencyMean + frequencyAlpha*(xp(2)-frequencyMean) ;
    amplitudesCosPredicted = xp(2+(1:nHarmonics))*amplitudeAlpha;
    amplitudesSinPredicted = xp(2+nHarmonics+(1:nHarmonics))*amplitudeAlpha;
    k = (1:nHarmonics).';
    H = [sum(-amplitudesCosPredicted.*k.*sin(k*phasePredicted) + amplitudesSinPredicted.*k.*cos(k*phasePredicted)) ...
         0 ...
         cos(k*phasePredicted).' ...
         sin(k*phasePredicted).'];
    Re = H*Pp*H' + R;                                         % Measurement error convariance matrix
    Kf = Pp*H'*inv(Re);                                       % Kalman gain
    yp = sum(amplitudesCosPredicted.*cos(k*phasePredicted) +...
             amplitudesSinPredicted.*sin(k*phasePredicted));  % Measurement prediction
    ep = y(n)-yp;                                             % Measurement prediction error
    xf = xp + Kf*(y(n)-yp);
    Pf = Pp - Kf*H*Pp;                                        % Filtered state error covariance matrix
    Pp = F*Pf*F' + Q;                                         % Predicted state error covariance matrix
    xp = F*xf;
    xp(1) = mod(xp(1),2*pi);                                  % Modulus operator to reduce round-off error    
    xp(2) = frequencyMean + frequencyAlpha*(xf(2)-frequencyMean);

    %--------------------------------------------------------------------------
    % Store Estimates
    %--------------------------------------------------------------------------
    phaseFiltered         = xf(1);
    frequencyFiltered     = xf(2);
    amplitudesCosFiltered = xf(2+(1:nHarmonics)); 
    amplitudesSinFiltered = xf(2+nHarmonics+(1:nHarmonics));
    
    Estimates.Filtered.phase        (n)   = phaseFiltered;
    Estimates.Filtered.frequency    (n)   = frequencyFiltered;
    Estimates.Filtered.amplitudesCos(n,:) = amplitudesCosFiltered.';
    Estimates.Filtered.amplitudesSin(n,:) = amplitudesSinFiltered.';

    Estimates.Predicted.y            (n)   = yp;
    Estimates.Predicted.error        (n)   = ep;
    Estimates.Predicted.phase        (n)   = phasePredicted;
    Estimates.Predicted.frequency    (n)   = frequencyPredicted;
    Estimates.Predicted.amplitudesCos(n,:) = amplitudesCosPredicted.';
    Estimates.Predicted.amplitudesSin(n,:) = amplitudesSinPredicted.';  
    
    if smoother
        %----------------------------------------------------------------
        % Store Variables for Smoothing
        %----------------------------------------------------------------    
        Estimates.Smoothed.jacobianFunction(n).state       = F;
        Estimates.Smoothed.jacobianFunction(n).measurement = H;
        Estimates.Smoothed.ErrorCovariance(n).measurement  = Re;
        Estimates.Smoothed.ErrorCovariance(n).state        = Pp;
    end
end

%==========================================================
% EKS recursions (Bryson-Frazier)
%==========================================================
if smoother
    %------------------------------------------------------
    % Variable Initialization
    %------------------------------------------------------
    l = zeros(nStateVariables,1);                          % Adjoint variable used for smoothing

    for n=nSamples:-1:2,
        %--------------------------------------------------
        % Recall Variables of Interest
        %--------------------------------------------------
        F   = Estimates.Smoothed.jacobianFunction(n).state;
        H   = Estimates.Smoothed.jacobianFunction(n).measurement;
        Re  = Estimates.Smoothed.ErrorCovariance(n).measurement;
        e   = Estimates.Predicted.error(n);
        Pp  = Estimates.Smoothed.ErrorCovariance(n-1).state;

        %--------------------------------------------------
        % Bryson-Frazier Smoothing Equations in terms of
        % predicted quantities
        %--------------------------------------------------
        xp  = [Estimates.Predicted.phase(n);           ...
               Estimates.Predicted.frequency(n);       ...               
               Estimates.Predicted.amplitudesCos(n,:)';...
               Estimates.Predicted.amplitudesSin(n,:)' ...
               ];        
        Rei = inv(Re);                                     % Error (innovation) covariance matrix
        Kp = (F*Pp*H')*Rei;                                % Kalman predicted gain
        Fp = F - Kp*H;                                     % Predicted state error transition matrix
        l  = Fp'*l + H'*Rei*e;                             % Adjoint variable
        xs = xp + Pp*l;                                    % Smoothed state estimate
        
        %--------------------------------------------------
        % Bryson-Frazier Smoothing Equations in terms of
        % filtered quantities
        %--------------------------------------------------
%         xf  = [Estimates.Filtered.phase(n);           ...
%                Estimates.Filtered.frequency(n);       ...               
%                Estimates.Filtered.amplitudesCos(n,:)';...
%                Estimates.Filtered.amplitudesSin(n,:)' ...
%                ];
%         Rei = inv(Re);                                     % Error (innovation) covariance matrix
%         Kp = (F*Pp*H')*Rei;                                % Kalman predicted gain
%         Fp = F - Kp*H;                                     % Predicted state error transition matrix
%         xs = xf + Pp*Fp'*l;                                % Smoothed state estimate
%         l  = Fp'*l + H'*Rei*e;                             % Adjoint variable
%         
        %--------------------------------------------------
        % Store smoothed estimates
        %--------------------------------------------------        
        phaseSmoothed         = xs(1);     
        frequencySmoothed     = xs(2);
        amplitudesCosSmoothed = xs(2+(1:nHarmonics));
        amplitudesSinSmoothed = xs(2+(nHarmonics+1:2*nHarmonics));
        
        ys = sum(amplitudesCosSmoothed.*cos(k*phaseSmoothed)+...
            amplitudesSinSmoothed.*sin(k*phaseSmoothed));  % Measurement prediction
        es = y(n)-ys;                                      % Measurement prediction error
   
        Estimates.Smoothed.y            (n)   = ys;
        Estimates.Smoothed.error        (n)   = es;
        Estimates.Smoothed.phase        (n)   = phaseSmoothed;        
        Estimates.Smoothed.frequency    (n)   = frequencySmoothed;
        Estimates.Smoothed.amplitudesCos(n,:) = amplitudesCosSmoothed;
        Estimates.Smoothed.amplitudesSin(n,:) = amplitudesSinSmoothed;
    end
    % In the absence of true smoothed estimates for the first sample
    % we will put the filtered estimate here instead of leaving these
    % as NaN or forcing them to be zero.
    Estimates.Smoothed.y            (1)   = Estimates.Filtered.y            (1);
    Estimates.Smoothed.error        (1)   = Estimates.Filtered.error        (1);
    Estimates.Smoothed.phase        (1)   = Estimates.Filtered.phase        (1);
    Estimates.Smoothed.frequency    (1)   = Estimates.Filtered.frequency    (1);
    Estimates.Smoothed.amplitudesCos(1,:) = Estimates.Filtered.amplitudesCos(1,:);
    Estimates.Smoothed.amplitudesSin(1,:) = Estimates.Filtered.amplitudesSin(1,:);
end


