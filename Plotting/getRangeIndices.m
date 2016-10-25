function [ rangeArray ] = getRangeIndices(sortedValues, rangeValues)
%GETRANGEINDICES Returns indices of values between min and max in sV
%   Assumes input array 'sortedValues' is sorted in descending order

maxValue = max(rangeValues);
minValue = min(rangeValues);
% Find the indices corresponding to the high-value and low-value limits
hiIndex = find(diff(sortedValues >= maxValue));
loIndex = find(diff(sortedValues >= minValue)); % diff returns arr len(n-1)

% If all values below highest value, will be no diff
if(isempty(hiIndex))
    hiIndex = 1;
end
% If all values above lowest value, will be no diff
if(isempty(loIndex))
    loIndex = length(sortedValues);
end

rangeArray = hiIndex:loIndex;

end

