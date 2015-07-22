function [model ...
          spontaneous_rate ...
          errors] = ConstrainedFit(input_data, ...
                                   target_data, ...
                                   spontaneous_rate, ...
                                   model, ...
                                   num_epochs, ...
                                   learning_rate, ...
                                   constrain, ...
                                   fix_idx)
%
%function [model
%          spontaneous_rate
%          errors] = ConstrainedFit(input_data,
%                                   target_data,
%                                   spontaneous_rate,
%                                   model,
%                                   num_epochs,
%                                   learning_rate,
%                                   constrain)
%
%   INPUT ARGUMENTS
%   input_data              The (N X M) input data matrix with N input rows
%                           of M parameters each
%   target_data             The (N X 1) target data matrix with N input rows
%   spontaneous_rate        These statements are equivalent:
%                              The spontaneous firing rate of the neuron
%                              The DC offset of the model output
%                              The input bias for the adaline neural
%                               network model
%                           Default: 0
%   model                   The M parameter linear model being fit
%                           Default = zeros(1,M)
%   num_epochs              The number of times to iterate over the entire
%                           set of training data.
%                           Default: 1000
%   learning_rate           The "gain" of the learning algorithm. Can be a
%                           single number for a constant learning rate, or
%                           can be presented as a (1 X L) vector for
%                           simulated annealing. In this case, num_epochs
%                           must also be a (1 X L) vector.
%                           Default: 0.001
%   constrain               If set to TRUE, constrains the criterion
%                           function such that:
%                           model error = sum((target - max(output,0)).^2)
%                           otherwise, it is an unconstrained linear
%                           optimization problem where:
%                           model error = sum((target - output).^2)
%                           Default: true
%
%   OUTPUT ARGUMENTS
%   model                   The trained model
%   spontaneous_rate        The DC offset of the model output.
%   errors                  The mean squared error of the model on the
%                           training data for each epoch.

%By default, don't penalize errors resulting from an estimated negative
%spike rate. 
if exist('constrain','var')
    if isempty(constrain)
        constrain = 1;
    end
else
    constrain = 1;
end

if ~exist('fix_idx','var')
    fix_idx = [];
end

total_rows = 0;
for data_cell_num = 1:length(input_data)
    input = input_data{data_cell_num};
    total_rows = total_rows + size(input,1);
end
model_order = size(input,2);
model_order = model_order + 1; %DC offset/Bias term

%Convert learning rate to [num_epochs1 num_epochs2 ... ; lr1 lr2 ...]
%where num_epochs* indicates how many epochs to use lr* in.
if length(learning_rate) == 1
    learning_rate = [num_epochs ; learning_rate];
end
if size(learning_rate,1) ~= 2
    error('Learning rate incorrectly specified');
end
learning_rate(1,:) = cumsum(learning_rate(1,:));
learning_rate(2,:) = learning_rate(2,:) / total_rows;
num_epochs = learning_rate(1,end);

model = model(:);
model = [model ; spontaneous_rate]; %DC offset/Bias term
errors = zeros(1,num_epochs);
if num_epochs > 0
    display('Performing Gradient Descent Optimization...');
    for epoch_num = 1:num_epochs
        if rem(epoch_num,100) == 0
            display([int2str(epoch_num) ' / ' int2str(num_epochs)]);
        end
        gradient = zeros(1,model_order);
        mse = 0;
        for data_cell_num = 1:length(input_data)
            input = input_data{data_cell_num};
            input = [input ones(size(input,1),1)];
            output = input * model;
            if constrain
                output(output<0) = 0;
            end
            model_error = target_data{data_cell_num} - output; 
            mse = mse + sum(model_error.^2);
            gradient = gradient + model_error'*input;
        end
        errors(epoch_num) = mse;
        [tmp lr_idx] = find(epoch_num <= learning_rate(1,:));
        lr = learning_rate(2,lr_idx(1));
        gradient(fix_idx) = 0;
        model = model + (lr*gradient)';
    end
end

spontaneous_rate = model(end);
model = model(1:end-1);
