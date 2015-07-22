function H = MatrixInformationTransfer(matrix)

% matrix = eye(13)
% matrix(1,1) = 2
% 
% matrix = ones(5,5)*0.1
% matrix = zeros(2,2);
% matrix(1,2) = 1;
% matrix(2,1) = 1;

[nRows nColumns] = size(matrix);

Ntotal = sum(sum(matrix));

H = 0;
for iRow=1:nRows
    for iColumn=1:nColumns
        if matrix(iRow,iColumn) ~= 0

            columnSum = 0;
            for jRow = 1:nRows
                columnSum = columnSum + matrix(jRow,iColumn);
            end
            rowSum = 0;
            for jColumn = 1:nColumns
                rowSum = rowSum + matrix(iRow,jColumn);
            end

            H = H + ...
                matrix(iRow,iColumn) ...
                * (log2(matrix(iRow,iColumn)) - ...
                   log2(columnSum) - ...
                   log2(rowSum) + ...
                   log2(Ntotal));
        end
    end
end

H = H/Ntotal;