function array_form = Cell2Array(cell_form)
%
%function array_form = Cell2Array(cell_form)
%
%   INPUT ARGUMENTS
%   cell_form           Data in Bat2Matlab's cell format
%
%   OUTPUT ARGUMENTS
%   array_form          Data in Matlab's array format
%
%Cell2Array converts data in Bat2Matlab's cell format 
%into Matlab's standard array format

array_form = [];
for data_cell_num = 1:length(cell_form)
    array_form = [array_form ; sparse(cell_form{data_cell_num})];
end