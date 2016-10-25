function [plotData, fitData] = ...
    plotTrends(sampleData, xVar, yVar, legendVars, varargin)
% PLOTTRENDS Plots trends in fitted impedance data
%---Inputs
% sampleData - vector of user-defined runs
% xVar - the name of the x-variable
% yVar - the index of the fit parameter to be used as the y-var
%   If yVar is given as a cell, e.g. {2, 'pO2'}, the second variable is
%   used as the timetrace variable
%   If yVar is given as array, e.g. [2 3], multiple fit params are plotted
% legendVars - string or cell array of fields to be used as the legend
% varargin - variable length options list
%  '(length|area|volume) normalized' - normalize by length, area, or vol
%  'flip colors' - reverse the order of colors
%  [-1 1 ...] - sets ascending (1) or descending (-1) legend var sort
%  'plot linear|loglog|semilogxy|arrhenius' sets plotting style
%  'plot tt linear|loglog|semilogxy' sets timetrace plotting style
%  'fit linear|loglog|semilogxy|arrhenius' fits and returns fit data in
%       given style
%  'average' - averages Y data
%  'errorbars' - plots error bars
%  'residuals' - makes a sub-plot of residuals
%  'write to file' - writes resulting data to tab-delimited text file


% Get rid of runs we don't want to include
for i = length(sampleData):-1:1
    if(sampleData(i).includePlot == 0)
        sampleData(i) = [];
    end
end
if(isempty(sampleData))
    error('Error: no data selected to include in plot.');
end

% Detect timetrace
if(iscell(yVar))
    ttVar = yVar{2};
    yVar = yVar{1};
else
    ttVar = '';
end

% Get legend indices
varNames = fieldnames(sampleData);
legendIndices = getStringIndex(legendVars, varNames);

% Get options
plotOptions = configureOptions(varargin, legendIndices);
if(strcmp(ttVar, ''))
    plotOptions.timetraceStyle = 0;
end

% Make the plot!
fh = figure('Color', 'White', 'Position', [100 100 900 600]);
[plotData, fitData] = plotParameters(sampleData, varNames, ...
    xVar, yVar, ttVar, legendVars, plotOptions);

% Write the data to a file
if(plotOptions.writeToFile)
    yVarName = sampleData(1).fitParamNames{yVar};
    sampleName = sampleData(1).sample;
    writeDataToFile(plotData, fitData, xVar, yVarName, sampleName, legendVars);
end

end

%% Set values for all the relevant options
function [opts] = configureOptions(args, lIndex)

% Set return values to defaults
norm = 0; fits = 0; style = 0; tt = 1; colorFlip = 0; bars = 0;
res = 0; avg = 0; lineStyle = '-'; order = lIndex; write = 0;
scatter_points = 0; cmap = 'viridis'; differentShapes = 0;

% Turn on options as the user requests
for i = 1:length(args)
    argVar = args{i};
    if(isnumeric(argVar))
        order = argVar .* order;
    else
        switch(argVar)
            case{'fit linear'}
                fits = 1;
            case{'fit semilogx'}
                fits = 2; style = 2;
            case{'fit semilogy'}
                fits = 3; style = 3;
            case{'fit loglog'}
                fits = 4; style = 4;
            case{'fit arrhenius'}
                fits = 5; style = 5;
            case{'plot linear'}
                style = 1;
            case{'plot semilogx'}
                style = 2;
            case{'plot semilogy'}
                style = 3;
            case{'plot loglog'}
                style = 4;
            case{'plot arrhenius'}
                style = 5;
            case{'plot tt linear'}
                tt = 1;
            case{'plot tt semilogx'}
                tt = 2;
            case{'plot tt semilogy'}
                tt = 3;
            case{'plot tt loglog'}
                tt = 4;
            case {'length norm', 'length normalized', 'length'}
                norm = 1;
            case {'area norm', 'area normalized', 'area'}
                norm = 2;
            case {'volume norm', 'volume normalized', 'volume'}
                norm = 3;
            case {'flip', 'flip colors', 'reverse', 'reverse colors'}
                colorFlip = 1;
            case {'bars', 'errorbars','error bars','errors'}
                bars = 1;
            case {'residuals'}
                res = 1;
            case {'avg', 'average'}
                avg = 1;
            case 'no lines'
                lineStyle = 'none';
            case 'write to file'
                write = 1;
            case 'scatter points'
                scatter_points = 1;
            case {'magma', 'plasma', 'inferno'}
                cmap = argVar;
            case {'change shapes', 'different shapes', 'shapes'}
                differentShapes = 1;
            otherwise
                error(['Option ' argVar ' unrecognized']);
        end
    end
end

% Create return struct with all options as fields
opts = struct('normalization', norm, ...
    'fitStyle', fits, ...
    'plotStyle', style, ...
    'timetraceStyle', tt, ...
    'colorFlip', colorFlip, ...
    'errorbars', bars, ...
    'residuals', res, ...
    'averaging', avg, ...
    'sortOrder', order, ...
    'lineStyle', lineStyle, ...
    'writeToFile', write, ...
    'scatterPoints', scatter_points, ...
    'cmap', cmap, ...
    'diffShapes', differentShapes);
end

%% Plot everything based on the input and the configuration
function [plotData, fitData] = ...
    plotParameters(structData, varNames, xVar, yVar, ttVar, lVars, opts)

% Transform the data from a struct to a cell array and sort by legend vars
xIndex = getStringIndex(xVar, varNames);
if(isempty(xIndex))
    error(['Couldn''t find x-variable ' xVar '.']);
end
[cellData, firsts, lasts, legendEntries] = ...
    getSortedStructData(structData, opts.sortOrder, xIndex);

% Extract the x- and y- data
[xData, yData, yErrData, yAxisLabel] = ...
    collectData(cellData, varNames, xIndex, yVar, ...
    opts.normalization, opts.plotStyle);

% Average points at same value of x together if averaging selected
if(opts.averaging)
    [xData, yData, yErrData, firsts, lasts] = ...
        averageData(xData, yData, yErrData, firsts, lasts);
end

% Collect information about the area of each point for the scatter plot
if(opts.scatterPoints)
    areaIndex = getStringIndex('area', varNames);
    areas = cell2mat(cellData(:, areaIndex));
end
shapes = {'^', 's', 'p', 'h'};

% Set up the figure, collect axis handle and colors array
[ah, colors] = setupFigure(length(firsts), opts.colorFlip, ...
    opts.residuals, opts.timetraceStyle, opts.cmap);

% Extract and plot timetrace data
if(opts.timetraceStyle)
    [tData, ttYData] = collectTimetraceData(cellData, varNames, 'hour', ttVar);
    ttPlot = plot(ah(2), tData, ttYData, 'k', 'LineWidth', 2);
    ylabel(ah(2), getLabelWithUnits(ttVar));
end

% Initialize saved slopes and intercepts
if(opts.fitStyle)
    slopes = [];
    intercepts = [];
end
plotData = {};
fitData = {};
hold on;
[nLegendRows, nLegendCols] = size(legendEntries);
for i = 1:length(firsts)
    color = colors(i, :); % Current color
    faceColor = color;
    if(nLegendCols > 1 && legendEntries(i, 2) == 0)
        faceColor = 'none';
    end
    shape = 'o';
    if(opts.diffShapes)
        shape = shapes{mod(i-1, length(shapes)) + 1};
    end
    
    dataName = mat2str(legendEntries(i, :));
    dataRange = firsts(i):lasts(i); % Current datarange
    plotData(i).xData = xData(dataRange);
    plotData(i).yData = yData(dataRange);
    plotData(i).yErrData = yErrData(dataRange);
    plotData(i).legend = dataName;
    if(opts.errorbars) % Error-bar plot
        lowerBar = min(yErrData(dataRange), yData(dataRange)*0.99999);
        upperBar = yErrData(dataRange);
        barPlot = errorbar(ah(1), xData(dataRange), yData(dataRange), ...
            lowerBar, upperBar, 'MarkerFaceColor', faceColor, ...
            'MarkerEdgeColor', color / 1.1, 'Color', color, 'Marker', shape, ...
            'LineStyle', opts.lineStyle, 'MarkerSize', 10);
    elseif(opts.scatterPoints)
        rawPlot = scatter(ah(1), xData(dataRange), yData(dataRange), ...
            areas(dataRange)*3e5, color, 'MarkerFaceColor', faceColor, ...
            'MarkerEdgeColor', color / 1.1, 'DisplayName', dataName, ...
            'Marker', shape);
    else % Regular plot
        rawPlot = plot(ah(1), xData(dataRange), yData(dataRange), ...
            'Color', color, 'MarkerFaceColor', faceColor, 'MarkerEdgeColor', ...
            color / 1.1, 'DisplayName', dataName, ...
            'Marker', shape, 'LineStyle', opts.lineStyle, 'MarkerSize', 10);
    end
    if(opts.fitStyle) % Plot fit data
        % Fit the plotted data
        [xFitData, yFitData, yRes, slope, intercept] = ...
            getFit(xData(dataRange), yData(dataRange), yErrData(dataRange),...
                opts.fitStyle);
        fitData(i).xData = xFitData;
        fitData(i).yData = yFitData;
        fitData(i).legend = dataName;
        fitData(i).slope = slope;
        fitData(i).intercept = intercept;
        fitPlot = plot(ah(1), xFitData, yFitData, ...
            'LineStyle', '-', 'Marker', 'none', 'Color', color, ...
            'LineWidth', 1);
        % Remove fit from legend
        set(get(get(fitPlot, 'Annotation'), 'LegendInformation'), ...
            'IconDisplayStyle', 'off');
        if(opts.residuals)
            resPlot = plot(ah(2), xData(dataRange), yRes, ...
                'LineStyle', opts.lineStyle, 'Marker', 'x', 'Color', color, ...
                'MarkerSize', 10, 'LineWidth', 3);
        end
    end
end
hold off;

% Finish figure by creating legend
lh = finishFigure(ah, opts.plotStyle, opts.residuals, ...
    yAxisLabel, xVar, legendEntries, lVars);

% If there's no fit, return an empty struct, otherwise print slopes
if(~opts.fitStyle)
    fitData = struct([]);
else
    lPos = get(lh, 'Position');
    displayText = {num2str(vertcat(fitData(:).slope), '%0.2g')};
    th = text(plotData(1).xData(1), plotData(1).yData(1), displayText);
    set(th, 'FontName', 'Arial', 'FontSize', getFontSize('slope'), ...
        'Units', 'Normalized', 'HorizontalAlignment', 'Center')
end

end

%% Collect the data from the cell into the correct matrix
function [xData, yData, yErrors, yAxisLabel] = ...
    collectData(cellData, varNames, xIndex, paramIndex, norm, style)

% Get indices for normalization and plotting parameters
normIndices = getStringIndex({'length', 'area', 'volume'}, varNames);
yIndices = getStringIndex({'fitParamValues', 'fitParamErrors', ...
    'fitParamNames'}, varNames);

% Extract data from cell matrix
xData = cell2mat(cellData(:, xIndex));
fitParamData = cell2mat(cellData(:, yIndices(1)));
fitParamErrors = cell2mat(cellData(:, yIndices(2)));
fitParamNames = cellData{1, yIndices(3)};
normData = cell2mat(cellData(:, normIndices));

fitElem = fitParamNames{paramIndex};

% Extract data and set normalization power (inverse for capacitance)
if(strncmpi(fitElem, 'CPE', 3)) % 'CPE' element
    power = -1;
    label = 'Effective capacitance / F';
    Y = fitParamData(:, paramIndex);
    n = fitParamData(:, paramIndex + 1);
    R = fitParamData(:, paramIndex - 1);
    % Compute 'effective capacitance' for a CPE
    yData = Y .^ (1 ./ n) .* R .^ ( (1./n) - 1);
    Yerr = fitParamErrors(:, paramIndex);
    nerr = fitParamErrors(:, paramIndex + 1);
    Rerr = fitParamErrors(:, paramIndex - 1);
    % Add errors in quadrature for a CPE
    dCeffdY = yData ./ (Y .* n);
    dCeffdn = -yData .* (log(R) + log(Y)) ./ (n .^ 2);
    dCeffdR = (1./n - 1) .* yData ./ R;
    yErrors = sqrt((Yerr .* dCeffdY).^2 + (nerr .* dCeffdn).^2 + ...
        (Rerr .* dCeffdR).^2);
elseif(strncmpi(fitElem, 'QPE', 3)) % 'QPE' element
    power = -1;
    label = 'Effective capacitance / F';
    Q = fitParamData(:, paramIndex);
    n = fitParamData(:, paramIndex + 1);
    R = fitParamData(:, paramIndex - 1);
    % Compute 'effective capacitnace' for a QPE
    yData = Q .* R .^ ( (1./n) - 1);
    Qerr = fitParamErrors(:, paramIndex);
    nerr = fitParamErrors(:, paramIndex + 1);
    Rerr = fitParamErrors(:, paramIndex - 1);
    dCeffdQ = yData ./ Q;
    dCeffdn = -yData .* (log(R)) ./ (n .^ 2);
    dCeffdR = (1./n - 1) .* yData ./ R;
    yErrors = sqrt((Qerr .* dCeffdQ).^2 + (nerr .* dCeffdn).^2 + ...
        (Rerr .* dCeffdR).^2);
elseif(strncmpi(fitElem, 'C', 1)) % Capacitor
    power = -1;
    label = 'Capacitance / F';
    yData = fitParamData(:, paramIndex); % yData is just capacitance
    yErrors = fitParamErrors(:, paramIndex);
elseif(strncmpi(fitElem, 'R', 1)) % Resistor
    power = 1;
    label = 'Resistance / \Omega';
    yData = fitParamData(:, paramIndex); % yData is just resistance
    yErrors = fitParamErrors(:, paramIndex);
else
    power = 1;
    label = fitElem;
    yData = fitParamData(:, paramIndex);
    yErrors = 0.5 * yData;
end

% Nomralize data appropriately
if(norm > 0)
    yData = yData .* (normData(:, norm) .^ power);
    yErrors = yErrors .* (normData(:, norm) .^ power);
end

% Arrhenius plotting
if(style == 5)
    xData = 1000 ./ (xData + 273.15);
    yData = yData; % log(yData); %%%%%%%%%
end

% If normalized by area, append appropriate suffix
suffix = '';
if(norm == 1)
    if(power < 0)
        suffix = 'cm^{-1}';
    else
        suffix = 'cm';
    end
elseif(norm > 1)
    if(power > 0)
        suffix = ['cm^' num2str(norm)];
    else
        suffix = ['cm^{-' num2str(norm) '}'];
    end
end

yAxisLabel = [label, suffix];

end

%% Collects data for timetrace plot
function [tData, yData] = collectTimetraceData(cellData, names, tVar, yVar)

% Get the indices and sort by time
indices = getStringIndex({tVar, yVar}, names);
ttData = cell2mat(cellData(:, indices));
ttData = sortrows(ttData, 1);
tData = ttData(:, 1);
yData = ttData(:, 2);
[changes, ~] = find(diff(yData));
changes = unique(changes);
% At each change, insert a step change into the timetrace
for i = 1:length(changes)
    j = changes(i) + i - 1;
    tData = [tData(1:j); tData(j); tData(j+1:end)];
    yData = [yData(1:j); yData(j+1); yData(j+1:end)];
end
end

%% Average the y points at the same x in the same legend group
function [xAvg, yAvg, yAvgErr, newFirsts, newLasts, nAvg] = ...
    averageData(xData, yData, yErrData, firsts, lasts)

% Collect the original data
dataMat = [xData, yData, yErrData];
avgData = zeros(length(dataMat(:, 1)), 3);
newFirsts = ones(length(firsts) + 1, 1);
newLasts = ones(length(firsts), 1);
nAvg = ones(length(firsts));
for i = 1:length(firsts)
    % Sort within each legend grouping
    curRange = firsts(i):lasts(i);
    nAvg(i) = length(curRange);
    sortedData = sortrows(dataMat(curRange, :), 1);
    % Get unique X within legend grouping
    [xData, xFirsts, ~] = unique(sortedData(:, 1));
    if(length(xFirsts) > 1)
        xLasts = [xFirsts(2:end) - 1; length(curRange)];
    else
        xLasts = length(curRange);
    end
    % Average data within each unique X within each grouping
    yData = zeros(length(xData), 1);
    yErrData = zeros(length(xData), 1);
    for j = 1:length(xFirsts)
        yData(j) = mean(sortedData(xFirsts(j):xLasts(j), 2)); % mean of values
        nDataPts = length(xFirsts(j):xLasts(j));
        %yErrData(j) = sqrt(sum(sortedData(xFirsts(j):xLasts(j), 3).^2)) / nDataPts; %quadrature
        if(nDataPts == 1)
            yErrData(j) = sortedData(xFirsts(j):xLasts(j), 3);
        else
            yErrData(j) = std(sortedData(xFirsts(j):xLasts(j), 2)); %stdev of values
        end
    end
    
    % Update the firsts and lasts arrays
    newLasts(i) = newFirsts(i) + length(xData) - 1;
    avgData(newFirsts(i):newLasts(i), :) = [xData, yData, yErrData];
    newFirsts(i+1) = newLasts(i) + 1;
end

avgData = avgData(1:newLasts(end), :); % Trim the avgData array
xAvg = avgData(:, 1); % Output points x with averaged pts y & errors yErr
yAvg = avgData(:, 2);
yAvgErr = avgData(:, 3);
newFirsts = newFirsts(1:end-1); % Trim new 'firsts' array

end

%% Fit the data based on plotting style
function [xFitData, yFitData, yRes, slope, intrcpt] = ...
    getFit(xData, yData, yErrData, fitType)

normalizedErrors = yErrData ./ (yData);
% Choose the type of fit
switch(fitType)
    case 2 % semilogx
        xData = log10(xData);
    case 3 % semilogy (also arrh)
        yData = log10(yData);
        yErrData = log10(yErrData);
    case 4 % loglog
        xData = log(xData); 
        yData = log(yData); 
        yErrData = log(yErrData);
    case 5 % Arrhenius
        yData = log(yData); %see line 332
        yErrData = log(yErrData);
end

% fitParams = polyfit(xData, yData, 1);
%normalizedErrors = yErrData ./ yData;
%zeroIndices = find(normalizedErrors == 0);
weights = 1 ./ ((normalizedErrors ) .^ 2);
%weights(zeroIndices) = 0;
%[weights(zeroIndices)] = deal(max(weights));
fitObject = fit(xData, yData, 'poly1', 'Weights', weights);
slope = fitObject.p1;
intrcpt = fitObject.p2;
yRes = yData;
xFitData = unique(xData);

switch(fitType)
    case 1 % linear
        xFitData = [0; xFitData];
        yFitData = slope * xFitData + intrcpt;
        yRes = yData - (slope * xData + intrcpt);
    case 2 % semilogx
        yFitData = slope * xFitData + intrcpt;
        xFitData = 10.^(xFitData);
        yRes = yData - (slope * xData + intrcpt);
    case 3 % semilogy (also arrh)
        xFitData = [0; xFitData];
        yFitData = 10.^(slope * xFitData + intrcpt);
        yRes = 10.^yData - 10.^(slope * xData + intrcpt);
    case 4 % loglog
        yFitData = exp(slope * xFitData + intrcpt);
        xFitData = exp(xFitData);
        yRes = exp(yData) - exp(slope * xData + intrcpt);
    case 5 % Arrhenius
        yFitData = exp(slope * xFitData + intrcpt);
        yRes = exp(yData) - exp(slope * xData + intrcpt);
end

end

%% Set up the figure properties
function [handles, colors] = setupFigure(numColors, flip, res, tt, cm)

% Set colors
switch(cm)
    case 'viridis'
        colors = viridis(numColors);
    case 'magma'
        colors = magma(numColors);
    case 'inferno'
        colors = inferno(numColors);
    case 'plasma'
        colors = plasma(numColors);
end
if(flip)
    colors = flipud(colors);
end
for i = 1:length(colors(:, 1))
    if(norm(colors(i, :)) > 1)
        colors(i, :) = colors(i, :) / 1.0 ;
    end
end

if(numColors == 1)
    colors = [0 0 0];
end

% If residual plot selected, set up a sub-plot
if(res)
    handles(2) = subplot(5, 1, 4:5);
    set(handles(2), 'YScale', 'linear', 'FontName', 'Arial', ...
        'FontSize', getFontSize('axis'));
    hold on;
    handles(1) = subplot(5, 1, 1:3);
elseif(tt)
    handles(2) = subplot(4, 1, 1);
    set(handles(2), 'YAxisLocation', 'Left', 'XTickLabel', [], ...
        'FontSize', getFontSize('axis'));
    switch(tt)
        case 2
            set(handles(2), 'XScale', 'log');
        case 3
            set(handles(2), 'YScale', 'log');
        case 4
            set(handles(2), 'XScale', 'log', 'YScale', 'log');
    end
    hold on;
    handles(1) = subplot(4, 1, 2:4);
else
    handles = gca;
end

for k = 1:length(handles)
    set([handles(k); findall(gca, 'Type','text')], 'FontName', 'Arial');
    set([handles(k); findall(gca, 'Type','text')], 'FontSize', getFontSize('axis'));
end

end

%% Write the data to a file
function writeDataToFile(plotData, fitData, xName, yName, sampleName, legendNames)
% Determine newline symbol
if(ismac)
    newline = '\n';
else
    newline = '\r\n';
end

% Create filename without overwriting
prefix = [sampleName '_' yName '_vs_' xName];
if(ischar(legendNames))
    legendDesc = [ '[' legendNames ']' ];
else
    for i = 1:length(legendNames)-1
        legendNames{i} = [legendNames{i} ','];
    end
    legendDesc = [ '[' horzcat(legendNames{:}) ']' ];
end
plotFilename = [prefix '_plot-data.txt'];
copies = 0;
while(exist(plotFilename, 'file'))
    copies = copies + 1;
    plotFilename = [prefix '_plot-data_' num2str(copies) '.txt'];
end

fitFilename = [prefix '_fit-data.txt'];
copies = 0;
while(exist(fitFilename, 'file'))
    copies = copies + 1;
    fitFilename = [prefix '_fit-data_' num2str(copies) '.txt'];
end


% Open plot data filename
dataFile = fopen(plotFilename, 'w');
legendData = {plotData.legend};
numGroups = length(legendData);
maxLength = 0; % Length of longest data column

% Write the header lines
header = ''; legendHeader = ''; legendValues = '';
for i = 1:numGroups
    header = [header xName '\t' yName '\tError\t'];
    legendHeader = [legendHeader legendDesc '\t' legendDesc '\t' legendDesc '\t'];
    legendValues = [legendValues plotData(i).legend '\t' plotData(i).legend '\t' plotData(i).legend '\t'];
    len = length(plotData(i).xData);
    if(len > maxLength)
        maxLength = len;
    end
end
header = [header newline];
legendHeader = [legendHeader newline];
legendValues = [legendValues newline];
fprintf(dataFile, [header legendHeader legendValues]);

% Write the actual data
for i = 1:maxLength
    dataStr = '';
    for j = 1:numGroups
        if(i <= length(plotData(j).xData))
            dataStr = [dataStr num2str(plotData(j).xData(i)) '\t'...
                num2str(plotData(j).yData(i)) '\t' ...
                num2str(plotData(j).yErrData(i)) '\t'];
        else
            dataStr = [dataStr ' \t \t \t'];
        end
    end
    dataStr = [dataStr newline];
    fprintf(dataFile, dataStr);
end

fclose(dataFile);

% Now do that for the fits
if(~isempty(fitData))
    fitFile = fopen(fitFilename, 'w');
    legendData = {fitData.legend};
    numGroups = length(legendData);
    maxLength = 0;
    header = ''; legendHeader = ''; legendValues = '';
    for i = 1:numGroups
        header = [header xName '\t' yName '\t'];
        legendHeader = [legendHeader legendDesc '\t' legendDesc '\t'];
        legendValues = [legendValues fitData(i).legend '\t' fitData(i).legend '\t' ];
        len = length(fitData(i).xData);
        if(len > maxLength)
            maxLength = len;
        end
    end
    header = [header newline];
    legendHeader = [legendHeader newline];
    legendValues = [legendValues newline];
    fprintf(fitFile, [header legendHeader legendValues]);
    
    for i = 1:maxLength
        dataStr = '';
        for j = 1:numGroups
            if(i <= length(fitData(j).xData))
                dataStr = [dataStr num2str(fitData(j).xData(i)) '\t'...
                    num2str(fitData(j).yData(i)) '\t'];
            else
                dataStr = [dataStr ' \t \t'];
            end
        end
        dataStr = [dataStr newline];
        fprintf(fitFile, dataStr);
    end
    fclose(fitFile);
    
end
end

%% Finish up the figure to look nice
function [legendHandle] = finishFigure(axisHandles, style, res, yAxisLabel, ...
    xAxisLabel, legendEntries, lVars)

if(res)
    axes(axisHandles(2));
    zeroLine = refline(0, 0);
    set(zeroLine, 'Color', 'black');
    uistack(zeroLine, 'bottom');
    set(axisHandles(1), 'XTickLabel', {});
    xlabel(axisHandles(2), getLabelWithUnits(xAxisLabel));
end
linkaxes(axisHandles, 'x');
% Set main axis as the first if there are multiple
axes(axisHandles(1));
if(~res)
    xlabel(getLabelWithUnits(xAxisLabel));
end
ylabel(yAxisLabel);
xlim('auto');
ylim('auto');
% Finish the figure properties
switch(style)
    case 1
        set(axisHandles(1), 'XScale', 'linear', 'YScale', 'linear');
    case 2
        set(axisHandles(1), 'XScale', 'log', 'YScale', 'linear');
        xLimits = xlim;
        xlim([10^floor(log10(xLimits(1)/1.1)), 10^ceil(log10(xLimits(2)*1.1))]);
    case 3
        set(axisHandles(1), 'XScale', 'linear', 'YScale', 'log');
        yLimits = ylim;
        ylim([10^floor(log10(yLimits(1)/1.1)), 10^ceil(log10(yLimits(2)*1.1))]);
    case 4
        set(axisHandles(1), 'XScale', 'log', 'YScale', 'log')
        xLimits = xlim;
        yLimits = ylim;
        xlim([10^floor(log10(xLimits(1))), 10^ceil(log10(xLimits(2)))]);
        ylim([10^floor(log10(yLimits(1))), 10^ceil(log10(yLimits(2)))]);
    case 5 % Arrhenius
        set(axisHandles, 'XScale', 'linear', 'YScale', 'log');
        axisHandles(2) = makeTopAxis(axisHandles);
end

% Make the legend
entries = mat2str(legendEntries);
if(length(legendEntries) > 1)
    entries = strsplit(entries(2:end-1), ';');
    entries = strrep(entries, ' ', ', ');
end
legendHandle = legend(gca, entries, 'Location', 'EastOutside');
set(legendHandle, 'FontName', 'Arial', 'FontSize', getFontSize('legend'));
legendPos = get(legendHandle, 'Position');

legendTitle = annotation('textbox', ...
    [legendPos(1), legendPos(2)+legendPos(4), 0.2, 0.1], ...
    'String', getLabelWithUnits(lVars), 'FontName', 'Arial', ...
    'FontSize', getFontSize('legend title'), 'HorizontalAlignment', 'center', ...
    'LineStyle', 'none', 'VerticalAlignment', 'middle');

titlePos = get(legendTitle, 'Position');
xPos = (legendPos(1) + legendPos(3) / 2) - (titlePos(3) / 2);
yPos = (legendPos(2) + legendPos(4));
set(legendTitle, 'Position', [xPos, yPos, titlePos(3), titlePos(4)]);

%set(lt, 'String', getLabelWithUnits(lVars), 'FontSize', getFontSize('title'), ...
%    'Units', 'normalized','FontName','Arial');

set(gcf, 'PaperPositionMode', 'auto');
if(length(axisHandles) > 1)
    set(axisHandles(1), 'Units', 'normalized');
    set(axisHandles(2), 'Units', 'normalized');
    pos1 = get(axisHandles(1), 'Position');
    pos2 = get(axisHandles(2), 'Position');
    set(axisHandles(1), 'Position', [pos1(1), pos1(2), 0.65, pos1(4)]);
    set(axisHandles(2), 'Position', [pos2(1), pos2(2), 0.65, pos2(4)]);
else
    pos = get(axisHandles, 'Position');
    set(axisHandles, 'Position', [pos(1), pos(2), 0.65, pos(4)]);
    box on;
    set(axisHandles, 'LineWidth', 1.5, 'TickLength', [0.02 0.05]);
end

end

%% Connect the top arrhenius axis
function [ax2] = makeTopAxis(ax1)
% Set up a new axis at the same position as the old one
ax2 = axes('position', get(ax1, 'position'), ...
    'XAxisLocation', 'top', ...
    'Color', get(ax1, 'Color'), ...
    'YTickLabel', [], ...
    'YAxisLocation', 'right', ...
    'YScale', get(ax1, 'YScale'), ...
    'YTick', get(ax1, 'YTick'), ...
    'XScale', get(ax1, 'XScale'));
set(ax1, 'Box', 'off', 'Color', 'none');
% Make sure the two axes scale together
linkaxes([ax1, ax2],'xy');
hlink = linkprop([ax2, ax1], {'Position','YTick','XScale','YScale','YMinorTick'});
setappdata(ax2, 'Axis_Linkage', hlink);
set(gcf, 'Children', [ax1 ax2]);
set(ax2, 'LineWidth', 1.5, 'TickLength', [0.02 0.05]);
set(ax1, 'LineWidth', 1.5, 'TickLength', [0.02 0.05]);
axes(ax1);
% Define the temperature range and step size
temps = (800:-25:600)';
% Set the correct labels at the positions defined by the arrhenius val
topTickPositions = (1000 ./ (temps + 273.15));
topTickLabels = num2str(temps);
set(ax2, 'XTick', topTickPositions, 'XTickLabel', topTickLabels, ...
    'FontName', 'Arial', 'FontSize', getFontSize('axis'));
xlabel(ax1, '1000 / T (1000 K^{-1})');
xlabel(ax2, 'Temperature / ^\circC');
end

