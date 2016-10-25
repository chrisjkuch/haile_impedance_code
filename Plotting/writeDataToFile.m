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