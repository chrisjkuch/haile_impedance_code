function plotNyquistData(sampleData, legendVars, freqRange, varargin)
% PLOTNYQUISTDATA
% INPUTS
%   sampleData - a struct with fields corresponding to run information
%   legendVars - a list of variables to group by in the legend
%   freqRange - range of frequencies to plot. 'all' if all desired
%   varargin - a list of options
%       'length normalized' - normalize by length
%       'area normalized' - normalize by area
%       'volume normalized' - normalize by volume
%       'flip colors' - reverse the order of colors
%       [-1 1 ...] - sets ascending (1) or descending (-1) leg var sort
%       'bode' - makes a bode plot

% Get variable indices
varNames = fieldnames(sampleData);
legendIndices = getStringIndex(legendVars, varNames);

% Configure options
%[norm, plotFits, flipColors, sortOrder] = ...
opts = configureOptions(varargin, legendIndices);

if(strcmp(freqRange, 'all'))
    freqRange = [0, Inf];
elseif(~isnumeric(freqRange) || length(freqRange) ~= 2)
    error('freqRange must be a numeric array or the string ''all''.');
end

% Set up the figure
figHandle = setupFigure(opts);

% Do the plotting
plotArcs(sampleData, freqRange, opts);

set(gcf, 'PaperPositionMode', 'auto');
end

function [plotOptions] = configureOptions(args, lIndex)

% Set return values to defaults
norm = 0; fits = 0; flipColors = 0; order = lIndex; bodePlot = 0;

% Turn on options as the user requests
for i = 1:length(args)
    argVar = args{i};
    if(isnumeric(argVar))
        order = argVar .* order;
    else
        switch(argVar)
            case{'fits', 'plot fits'}
                fits = 1;
            case {'length norm', 'length normalized', 'length'}
                norm = 1;
            case {'area norm', 'area normalized', 'area'}
                norm = 2;
            case {'volume norm', 'volume normalized', 'volume'}
                norm = 3;
            case {'flip', 'flip colors', 'reverse', 'reverse colors'}
                flipColors = 1;
            case {'bode', 'plot bode'}
                bodePlot = 1;
            otherwise
                error(['Option ' argVar ' unrecognized']);
        end
    end
end

plotOptions = struct(...
    'plotFits', fits, ...
    'normalization', norm, ...
    'flipColors', flipColors, ...
    'sortOrder', order, ...
    'bodePlot', bodePlot);
end


function plotArcs(structData, freqs, opts)

% Remove runs we don't want to include
for i = length(structData):-1:1
    if(structData(i).includePlot == 0)
        structData(i) = [];
    end
end

if(isempty(structData))
    error('Error: no data selected to include in plot.');
end

% Get the sorted struct data
[cellData, firsts, lasts, legendEntries] = ...
    getSortedStructData(structData, opts.sortOrder, 0);

% Set colors
colors = colormap(viridis(length(firsts)));
if(opts.flipColors)
    colors = flipud(colors);
end
if(length(firsts) == 1)
    colors = [0, 0, 0];
end

% Get the cell array index of the impedance values
varNames = fieldnames(structData);
Zindices = getStringIndex({'Z', 'Zfit'}, varNames);
normIndices = getStringIndex({'length', 'area', 'volume'}, varNames);
if(opts.normalization > 0)
    normIndex = normIndices(opts.normalization);
else
    normIndex = 0;
end

% Plot arcs with appropriate legend grouping
% Plots each raw arc as dots, each fit with solid lines
hold all;
for i = 1:length(firsts)
    newEntry = 1; % Track if we're in a new legend group
    color = colors(i, :);
    
    for j = firsts(i):lasts(i)
        % Calculate normalization
        if(normIndex > 0)
            normFactor = cellData{j, normIndex};
        else
            normFactor = 1;
        end
        dataName = mat2str(legendEntries(i, :));
        xData = cellData{j, Zindices(1)}.re * normFactor;
        yData = -cellData{j, Zindices(1)}.im * normFactor;
        wData = cellData{j, Zindices(1)}.freq;
        zData = sqrt(xData.^2 + yData.^2);
        if(opts.bodePlot)
            xData = wData;
            yData = zData;
        end
        plottedData(i) = struct('xData', xData, 'yData', yData, ...
            'wData', wData, 'zData', zData);
        % Extract desired frequency range
        fRange = getRangeIndices(wData, freqs);
        
        rawArc = plot(xData(fRange), yData(fRange), ...
            'Color', color, 'Marker', 'o', 'LineStyle', 'none', ...
            'MarkerFaceColor', color, 'MarkerEdgeColor', color / 1.1, ...
            'MarkerSize', 6);
        set(rawArc, 'DisplayName', dataName);
        %div10Indices = find(mod(wData(fRange), 10) == 0);
        %rawArcDecades = plot(xData(div10Indices), yData(div10Indices), ...
        %    'Color', color, 'Marker', 's', 'LineStyle', 'none', ...
        %    'MarkerFaceColor', color / 1.1, 'MarkerEdgeColor', color / 1.1, ...
        %    'MarkerSize', 8);
        %set(get(get(rawArcDecades, 'Annotation'), 'LegendInformation'), ...
        %        'IconDisplayStyle', 'off');
        if(~newEntry)
            set(get(get(rawArc, 'Annotation'), 'LegendInformation'), ...
                'IconDisplayStyle', 'off');
        end
        % Plot fits
        if(opts.plotFits)
            xFitData = cellData{j, Zindices(2)}.re * normFactor;
            yFitData = -cellData{j, Zindices(2)}.im * normFactor;
            wFitData = cellData{j, Zindices(2)}.freq;
            zFitData = sqrt(xFitData.^2 + yFitData.^2);
            if(opts.bodePlot)
                xFitData = wFitData;
                yFitData = zFitData;
            end
            fRange = getRangeIndices(wFitData, freqs);
            
            fitArc = plot(xFitData(fRange), yFitData(fRange), ...
                'Color', color, 'LineWidth', 1, 'Marker', 'none');
            set(get(get(fitArc, 'Annotation'), 'LegendInformation'), ...
                'IconDisplayStyle', 'off');
            set(fitArc, 'DisplayName', dataName);
        end
        newEntry = 0;
    end
end
hold off

% Make the legend
entries = mat2str(legendEntries);
if(length(legendEntries) > 1)
    entries = strsplit(entries(2:end-1), ';');
    entries = strrep(entries, ' ', ', ');
end

legendHandle = legend(gca, entries, 'Location', 'EastOutside', ...
    'LineWidth', 1.5);
set(legendHandle, 'FontSize', getFontSize('legend'));
legendPos = get(legendHandle, 'Position');

legendTitle = annotation('textbox', ...
    [legendPos(1), legendPos(2)+legendPos(4), 0.4, 0.2], ...
    'String', getLabelWithUnits(varNames(abs(opts.sortOrder))), 'FontName', 'Arial', ...
    'FontSize', getFontSize('legend title'), 'HorizontalAlignment', 'center', ...
    'LineStyle', 'none', 'VerticalAlignment', 'middle', 'FitBoxToText', 'on');

titlePos = get(legendTitle, 'Position');
xPos = (legendPos(1) + legendPos(3) / 2) - (titlePos(3) / 2);
yPos = (legendPos(2) + legendPos(4));
set(legendTitle, 'Position', [xPos, yPos, titlePos(3), titlePos(4)]);

% Set X and Y limits  0-max
xLimits = xlim;
yLimits = ylim;

box on;
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1), pos(2), 0.62, pos(4)]);
end

% Set up the figure to look nice
function [figHandle] = setupFigure(opts)
figHandle = figure('Color', 'white', 'Position', [100, 100, 900, 600]);
set(gca, 'FontName', 'Arial',...
    'FontSize', getFontSize('axis'), ...
    'LineWidth', 1.5, 'TickLength', [0.02 0.05]);
suffix = '';
norm = opts.normalization;
if(norm == 1)
    suffix = '-cm';
elseif(norm == 2)
    suffix = '-cm^2';
elseif(norm == 3)
    suffix = '-cm^3';
end
bodePlot = opts.bodePlot;
if(bodePlot)
    xlabel('Frequency / Hz');
    ylabel(['Impedance / \Omega' suffix]);
    set(gca, 'XScale', 'log', 'YScale', 'log');
else
    xlabel(['Re(Z) / \Omega' suffix]);
    ylabel(['-Im(Z) / \Omega' suffix]);
    set(gca, 'DataAspectRatio', [1 1 1]);
end
end
