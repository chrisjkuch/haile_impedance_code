function [ labels ] = getLabelWithUnits( varNames )
%GETLABELWITHUNITS Returns string of label with units
%   Put in a variable name, get a label with units out

labels = '';
if(iscell(varNames))
    for i = 1:length(varNames)
        labels{i} = getUnits(varNames{i});
    end
else
    labels = getUnits(varNames);
end

end

function [unit] = getUnits(name)
unit = '';
switch(name)
    case {'diam', 'diameter'}
        unit  = 'Diameter / \mum';
    case 'Tstage'
        unit = 'Stage Temperature / ^\circC';
    case 'Tfilm'
        unit = 'Film Temperature / ^\circC';
    case 'Tsample'
        unit = 'Sample temperature / ^\circC';
    case 'pO2'
        unit = 'pO_2 / atm';
    case 'hour'
        unit = 'Time / h';
    case 'Tgrowth'
        unit = 'Growth Temperature / ^\circC';
    case 'thickness'
        unit = 'Thickness / nm';
    case 'bias'
        unit = 'Bias / mV';
    case 'indexI'
        unit = 'Row index';
    case 'indexJ'
        unit = 'Column index';
    case 'seqNum'
        unit = 'Run #';
    case {'xA', 'Astoich'}
        unit = 'A-site stoichiometry';
    case 'repeat'
        unit = 'Repeat run';
    case 'humidity'
        unit = '% water';
    case 'length'
        unit = 'Three-phase boundary length / cm';
    case 'area'
        unit = 'Two-phase boundary area / cm^2';
    case 'distance'
        unit = 'Distance from paste / mm';
    case 'thermovoltage'
        unit = 'Thermovoltage / mV';
    case 'xPos'
        unit = 'Distance from left / mm';
    case 'yPos'
        unit = 'Distance from top / mm';
    case 'cycle'
        unit = 'Cycle / #';
    case 'compositionLSCF'
        unit = 'x in La_{0.6}Sr_{0.4}Co_{1-x}Fe_{x}O_{3-\delta}';
    case 'compositionPBSCF'
        unit = 'x in PrBa_{0.5}Sr_{0.5}Co_{2-x}Fe_{x}O_{5+\delta}';
    case 'dotNum'
        unit = 'Dot number';
    case 'orientation'
        unit = 'hkl';
    otherwise
        unit = name;
        %error(['Error: no unit specified for variable named ' name '.']);
end

end