function numList = evalStrList(str)
%returns a list of numbers from a given string

numList = [];
while ~isempty(str)
    try
        [num str] = strtok(str);
        num = eval(num);
        numList = [numList num];
    catch e
        disp(['! invalid trace number input, ' num])
        continue
    end
end