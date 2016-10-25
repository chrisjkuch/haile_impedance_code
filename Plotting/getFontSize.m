function [ size ] = getFontSize( varargin )
%GETFONTSIZE Returns font size to use

if(ismac)
    size = 30;
else
    if(length(varargin) == 1)
        switch varargin{1}
            case 'legend'
                size = 16;
            case 'title'
                size = 20;
            case 'axis'
                size = 20;
            case 'slope'
                size = 16;
            case 'legend title'
                size = 16;
        end
    else
        size = 20;
    end
end    

end

