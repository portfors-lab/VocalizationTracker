function [varargout] = VarianceSeries(varargin)
%VarianceSeries: Time-series of the variance of signal intervals.
%
%   [y] = VarianceSeries(x,fs,wl,pf);
%
%   x    Input Signal.
%   fs   Sampling frequency (Hz).
%   wl   Window length,default = 100 samples
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y   A vector of standard deviation values.
%
%   The VarianceSeries function calculate the Variances of
%   succesive intervals of a time-series and returns a vector of 
%   these Variance values. Convolution with the fast Fourier 
%   transform is used to speed the calculation
%
%   Example: Calculate the variances of 2500 sample windows of a
%   electrocardiogram signal, and plot result.
%
%      load ECG;
%      VarianceSeries(ecg,fs,2500,1);
%
%   References
%
%   Version 0.01.01.21 DT
%   
%   See also CONV.

%=====================================================================
% Change log
%=====================================================================
%
%   Date         :  Description
%   ------------------------------------------------------------------
%   Oct 19, 2002 : Change from MATLAB conv function for convolution 
%                : to DT's FFTConvolution function
%   Oct 19, 2002 : Added row to column conversion rather than 
%                : generating an error
%                :
%=====================================================================


%=====================================================================
% Proccess Arguments
%=====================================================================


    if nargin < 2 | nargin > 4 | nargout > 1
        
        help('VarianceSeries')
        return
        
    end
    
    % extract input arguments
    x         = varargin{1};
    fs        = varargin{2};
    
    if nargin >= 3
        wl = varargin{3};
    else
        wl = 100;
    end
    
    if nargin == 4
        pf = varargin{4};
    else  
        if nargout == 0
            pf = 1;
        else
            pf = 0;
        end
    end
    
    % check if x is column vector
    x = x(:);
        
%=====================================================================
% Error Checking
%=====================================================================
    
    if length(x) == 0
        error('x is empty.');
        return
    elseif wl == 100 & length(x) < 100
        error('If wl is not specified, length(x) must be > 100');
        return
    elseif isa(x,'numeric') ~= 1  
        error('x must be numeric');
        return
   elseif length(fs)~=1
        error('fs must be scalar.');
        return
    elseif isa(fs,'numeric') ~= 1 
        error('fs must be numeric');
        return
    elseif fs <= 0
        error('fs must be positive.');
        return
    elseif floor(fs) ~= fs
        error('fs must be an integer.');
        return
   elseif length(wl)~=1
        error('wl must be scalar.');
        return
    elseif isa(wl,'numeric') ~= 1 
        error('wl must be numeric');
        return
    elseif wl <= 0
        error('wl must be positive.');
        return
    elseif floor(wl) ~= wl
        error('wl must be an integer.');
        return
    elseif length(pf) ~= 1  
        error('pf must be scalar.');
        return
    elseif isa(pf, 'numeric') ~= 1
        error('pf must be numeric');
        return
    elseif pf ~= 0 & pf ~= 1
        error('pf must be a Boolean Value.');
        return
    end

%=====================================================================
% Variables
%=====================================================================


    % hopefully these will be allocated on stack
    % no dynamic allocation
    lenX          = length(x);               % signal length
    xSquared      = zeros(1,lenX);           % signal squared
    xSquaredBar   = zeros(1,wl+lenX-1);      % square of means
    xBarSquared   = zeros(1,wl+lenX-1);      % mean of squares
    window        = ones(wl,1)/wl;           % window to calc means
    
%=====================================================================
% variance
%=====================================================================

    % square signal
    xSquared = x.^2;
    
    % calc xBarSquared
    xBarSquared = FFTConvolution(x,window).^2;
    
    % Calculate xSquaredbar
    xSquaredBar = FFTConvolution(xSquared,window);
    
    y =  (xSquaredBar(wl:(lenX))-xBarSquared(wl:(lenX)))*(wl/(wl-1));
          
%=====================================================================
% Plotting
%=====================================================================

    if pf
        
        % lots of cryptic plot formatting commands
        figH  = figure;
        t     = (0:length(y)-1)/fs+wl/(fs*2);
        plotH = plot( t, y, 'r' );
        axH   = get(figH,'CurrentAxes');
        titH  = get(axH,'Title');
        ylabH = get(axH,'YLabel');
        xlabH = get(axH,'XLabel');
        titstr = sprintf('Variances of Intervals of length %3.3f s',wl/fs);
        set(titH,'String',titstr);
        set(ylabH,'String','Variance');
        set(xlabH,'String','Time (s)');
        
        FigureSet(1,6.5,4.017);

    end

%=====================================================================
% Output
%=====================================================================

  if nargout == 1
      varargout{1} = y;
  end
  
%=====================================================================