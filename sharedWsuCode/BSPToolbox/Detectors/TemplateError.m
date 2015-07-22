function [varargout] = TemplateError(varargin)
%TemplateError: Calculate template error.
%
%   [y] = TemplateError(x,fs,template,weight,pf);
%
%   x          Input Signal.
%   fs         Sampling frequency (Hz).
%   template   Template.
%   weight     A weighted window which is the same length as the 
%              template, default = ones(length(template),1).
%   pf         Plot flag: 0=none (default), 1=screen.
%
%   y          The weighted mean squared error between the template 
%              and the input signal x
%
%   The TemplateError function is a feature detection 
%   tool. The template represents the feature one wishes to detect
%   the presence of in a given signal. The template is compared
%   with the signal along sequential intervals and the difference
%   between the signal and the template is calculated, squared, 
%   multiplied by a weighted window squared and then summed for
%   each comparable interval in the signal. Convolution using the
%   fast Fourier transform is used to speed up the calculation.
%
%   Once the time-series of error values has been calculated, it can
%   be used to measure the degree that the signal matches the feature
%   corresponding to the template at different points. The output
%   signal is indexed such that the first element corresponds to the
%   error between the weighted template and the beginning segment of 
%   the input signal of the same length as the template. The length of
%   the output signal is equal to the length of the input signal minus 
%   the length of the template plus one and the edges are ignored. The
%   error is only calculated between the template and segments of the
%   signal which are the same length as the template, the edges are
%   ignored.
%
%   Example: Use TemplateError to dectect ECG spikes.
%
%      load NoisyECG;
%      ri       = ceil(rand(30,1)*(length(ei)-1));
%      int      = ceil((median(ei(2:length(ei)) - ...
%                           ei(1:length(ei)-1)))/2);
%      temp     = zeros(2*int+1,1);
%      for k    = 1:length(temp)
%       off     = k-int-1; 
%       temp(k) = mean(ecg(ei(ri)+off)); 
%      end;
%      wind     = (0.995)*ones(length(temp),1);
%      err      = TemplateError(ecg,fs,temp,wind,1);
%      si       = find(err <= min(error)*5)+int;
%      figure;
%      plot((0:length(ecg)-1)/fs,ecg,'g',(si-1)/fs,ecg(si),'.r');
%
%   References
%
%   Version 0.05.00.21 DT
%
%   See also Detectors and CONV.

%=====================================================================
% Proccess Arguments
%=====================================================================
    if nargin < 3 | nargin > 5 | nargout > 1        
        help('TemplateError')
        return;
        end
    
    % extract input arguments
    x         = varargin{1};
    fs        = varargin{2};
    template  = varargin{3};
    
    if nargin >= 4
        weight = varargin{4};
    else
        weight = ones(length(template),1);
    end
    
    if nargin == 5 
        pf = varargin{5};
    else  
        if nargout == 0
            pf = 1;
        else
            pf = 0;
        end
    end
    
    % make sure vectors are column vectors
    x        = x(:);
    template = template(:);
    weight   = weight(:);
        
%=====================================================================
% Error Checking
%=====================================================================        

    sigSize    = size(x);
    tempSize   = size(template);
    weightSize = size(weight);
    
    if sigSize(1) == 0
        error('x is empty.');
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
    elseif tempSize(1) == 0
        error('template is empty.');
        return
    elseif isa(template,'numeric') ~= 1
        error('template must be numeric');
        return
    elseif weightSize(1) == 0
        error('weight is empty.');
        return
    elseif isa(weight,'numeric') ~= 1
        error('weight must be numeric');
        return
    elseif weightSize(1) ~= tempSize(1)
        error('weight and template must be the same length');
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
    len_t          = length(template);       % length of template
    len_x          = length(x);              % signal length   
    Constant_part  = zeros(1,len_t);         % Constant part of sum
    Product_part   = zeros(1,len_t+len_x-1); % Product part of sum
    Squared_part   = zeros(1,len_t+len_x-1); % Squared part of sum
    
%=====================================================================
% Calculate the error
%=====================================================================

    % square weight function
    weight = weight.^2;

    % square template sum with weight function 
    Constant_part = sum( template.^2.*weight);
  
    % reverse template and weight
    weight   = weight(len_t:-1:1);
    template = template(len_t:-1:1);
	
    % set template equal to weight * template
    template = weight.*template;
    
    % calc product part
    Product_part = 2*FFTConvolution(x,template);
    
    % Calculate squared part
    x = x.^2;
    Squared_part = FFTConvolution(x,weight);
    
    y =  (Squared_part(len_t:(len_x)) - Product_part(len_t:(len_x)) + ... 
          Constant_part)/len_t;
          
%=====================================================================
% Plotting
%=====================================================================

    if pf

        % lots of cryptic plot formatting commands
        figH  = figure;
        t     = (0:length(y)-1)/fs;
        plot( t, y, 'r' );
        axH   = get(figH,'CurrentAxes');
        titH  = get(axH,'Title');
        ylabH = get(axH,'YLabel');
        xlabH = get(axH,'XLabel');
        set(titH,'String','Weighted Mean Squared Error');
        set(ylabH,'String','Error');
        set(xlabH,'String','Time (s)');  
        

    end

%=====================================================================
% Output
%=====================================================================

  if nargout == 1
      varargout{1} = y;
  end
  
%=====================================================================