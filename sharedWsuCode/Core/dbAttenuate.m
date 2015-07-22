function attenuatedSignal = dBAttenuate(signal, dBAttenuation)
%
%function attenuatedSignal = dBAttenuate(signal, dBAttenuation)
%
%   INPUT ARGUMENTS
%   signal              The signal to be attenuated
%   dbAttenuation       The amount of attenuation to apply.
%
%   OUTPUT ARGUMENTS
%   attenuatedSignal    The attenuated signal
%
%dBAttenuate attenuates the power of the input signal by the
%specified amount. The definition of dB used for the calculation is:
%
%   dB = 20 log(A1/A2)
%
%where A1 is the peak amplitude of the attenuated signal and A2 is
%the peak amplitude of the input signal signal (which is equal
%to 1).

ratioMetricAttenuation = exp(-dBAttenuation * log(10) / 20);
attenuatedSignal = signal * ratioMetricAttenuation;