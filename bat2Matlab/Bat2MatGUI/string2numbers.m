function numList = string2numbers(str)

    strList = textscan(str, '%s');
    strList = strList{1}.'; %turn into row vector
    numList = [];
    for strnum = strList
        num = str2double(strnum{1});
        if isnan(num)
            disp(['!Test ' strnum{1} ' : invalid test number']);
        else
            numList = [numList num];
        end
    end
end