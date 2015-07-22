%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File Name   : FFTConvolution.m                          %
%  Author      : Daniel Gilbert                            %
%  Date        : July 22, 2002                             %
%  Purpose     : BSP                                       %
%  Description : Calculates the convolution of the two     %
%                signals given as arguements using the     %
%                fast Fourier Transform                    %
%                                                          %
%  Arguments                                               %
%                                                          %
%       x,y  : Signals                                     %
%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conv_sum = FFTConvolution(x,y)

    % calculate span
    span = length(x)+length(y)-1;
    pad  = 2^ceil(log2(span));

    % calculates FFTs
    X = fft(x,pad);
    Y = fft(y,pad);
    
    % F = XY
    F = X.*Y;

    % use inverse FFT 
    conv_sum = real(ifft(F));