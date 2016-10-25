function [ index ] = getStringIndex( searchStrings, names )
%GETSTRINGINDEX Returns the indices of searchStrings in names
%---Inputs
%   searchStrings - cell array of strings to be searched for in names
%   names - cell array of strings
%---Outputs
%   index - scalar array of indices of matched strings

% Get the column index of each legend variable
if(iscell(searchStrings))
    index = zeros(length(searchStrings), 1);
    for i = 1:length(index)
        cur_index = find(strcmp(names, searchStrings{i}));
        if(isempty(cur_index))
            errmsg = ['Could not find string ''' searchStrings{i} ...
                ''' in names array.'];
            error(errmsg);
        else
            index(i) = cur_index;
        end
    end
else
    index = find(strcmp(names, searchStrings));
end


end

