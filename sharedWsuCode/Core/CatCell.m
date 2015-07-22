function new_cell = CatCell(cell1, cell2)
%
%function new_cell = CatCell(cell1, cell2)
%
%   INPUT ARGUMENTS
%   cell1               The first cell in the concatenation
%   cell2               The second cell in the concatenation
%
%   OUTPUT ARGUMENTS
%   new_cell            The concatenated cell
%
%CatCell concatenates two data structures in the Bat2Matlab
%data cell format.

new_cell = [];
for data_cell_num = 1:length(cell1)
    new_cell{data_cell_num} = [cell1{data_cell_num} cell2{data_cell_num}];
end