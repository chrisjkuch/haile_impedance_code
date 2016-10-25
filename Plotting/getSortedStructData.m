function [ cellData, firsts, lasts, legendEntries ] = ...
    getSortedStructData( structData, sortBy, xIndex )
%GETSORTEDSTRUCTDATA Returns cell array of struct data sorted by field
%--- Inputs
%   structData - a struct
%   sortBy - a scalar array of the indicies of vars to sort by
%--- Outputs
%   cellData - sorted cell matrix with struct fields as columns
%   firsts - array of indices of first legend group entry
%   lasts - array of indices of last legend group entry

% Take the transpose of the cell matrix from struct2cell to get
% struct fields along columns
cellData = struct2cell(structData)';
if(xIndex == 0)
    cellData = sortrows(cellData, sortBy);
else
    cellData = sortrows(cellData, [sortBy; xIndex]);
end
legendData = cell2mat(cellData(:, abs(sortBy)));

% Find the legend groupings
[changedIndices, ~] = find(diff(legendData));
changedIndices = unique(changedIndices);
if(isempty(changedIndices))
    firsts = [1];
    lasts = [length(legendData(:, 1))];
else
    firsts = [1; changedIndices + 1];
    lasts = [changedIndices; length(legendData(:, 1))];
end
legendEntries = legendData(firsts, :);

end