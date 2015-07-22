function output = ProcessModel(model, ...
                                       input)
%
%function output = ProcessPoleZeroModel(model,
%                                       input)
%
%   INPUT ARGUMENTS
%   model                   The model to process
%   input                   Input data for the model
%
%   OUTPUT ARGUMENTS
%   output                  The output of the model

output = [];
for data_cell_num = 1:length(input)
    input_data = input{data_cell_num};
    output = [output ; input_data*model.parameters];
end
output = output + model.spontaneous_rate;
