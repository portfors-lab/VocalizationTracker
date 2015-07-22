function numList = evalStrList2(str)
%returns list of numbers from given string in a cell array

numList = {};
index = 1;
while ~isempty(str)
    try
        [num str] = strtok(str);
        num = eval(num);
        numList{index} = num;
    catch e
        disp(['! invalid trace number input, ' num])
        continue
    end
    index = index+1;
end