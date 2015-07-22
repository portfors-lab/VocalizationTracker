function fit_parms = FitData(input,output)

input_sequence_length = size(input,2);
output_length = size(output,2);
data = [input output];
data_corr = cov(data);

R = data_corr(1:input_sequence_length,1:input_sequence_length);
d = data_corr(1:input_sequence_length,input_sequence_length+1:input_sequence_length+output_length);

fit_parms = inv(R)*d;