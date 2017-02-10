function [medianFitParams] = ...
    customImpedanceFit(runs, circuitVersion, startingParams, freqRange, ...
    sampleDir, fitFolderName)
% CUSTOMIMPEDANCEFIT Fits impedance data to a user-defined circuit
% INPUTS
%   runs            vector of run structs
%   circuitVersion  string of version of circuit to use
%   startingParams  Starting points for fit parameters
%   freqRange       Two-element array indicating min and max of freq range
%                   Input 'all' for all frequencies
%   folderName      Full path to where sim files are written
%   name            Filename of summary to write fit summary file
% OUTPUTS
%   medianFitParams      - fitted parameters

%% Set up lsqnonlin parameters

% Set the fitting function options
opts = optimoptions('lsqnonlin', ...
    'TolX', 1E-15, ...
    'TolFun', 1E-15, ...
    'MaxIter', 1000, ...
    ...'PlotFcns', {@optimplotx, @optimplotstepsize}, ...
    ...'Display', 'off', ...
    'MaxFunEvals', 50000);

weighting = 'modulus'; % 'unit', 'proportional', 'modulus'
writeToFile = 1;
makePlots = 0;
verbose = 1;

if(strcmp(freqRange, 'all'))
    freqRange = [0, Inf];
end
if(ismac)
    slash = '/';
    newline = '\n';
else
    slash = '\';
    newline = '\r\n';
end
%% Choose the circuit
switch(circuitVersion)
    case 'R'
        pNames      = {'R'};
        lowerBounds = 1;
        upperBounds = Inf;
    case 'RL'
        pNames      = {'R', 'L'};
        lowerBounds = [0 0];
        upperBounds = [Inf Inf];
    case 'RRQ'
        pNames         = {'R0', 'R1',  'Y1', 'n1'};
        lowerBounds    = [   1,    1,     0,  0.0];
        upperBounds    = [ Inf,  Inf,   Inf,  1.0];
    case 'RRQRQ'
        pNames         = {'R0', 'R1',  'Y1', 'n1',  'R2',  'Y2', 'n2'};
        lowerBounds    = [   1e2,    1,     0,  0.0,     1,     0,  0.0];
        upperBounds    = [ Inf,  Inf,   Inf,  1.0,   Inf,   Inf,  1.0];
    case 'RRQRQRQ'
        pNames         = {'R0', 'R1',  'Y1', 'n1',  'R2',  'Y2', 'n2', 'R3', 'Y3', 'n3'};
        lowerBounds    = [   1e2,    1,    0,  0.0,      1,     0,  0.0,    1,    0,  0.0];
        upperBounds    = [ Inf,  Inf,  Inf,  1.0,    Inf,   Inf,  1.0,  Inf,  Inf,  1.0];
    case '6a'
        pNames         = {'Rion', 'Rion_s', 'Cion_s', 'Cchem', 'Ceon_p', 'R0'};
        lowerBounds    = [     0,        0,        0,       0,        0,    0];
        upperBounds    = [   Inf,      Inf,      Inf,     Inf,      Inf,  Inf];
    case '7b'
        pNames         = {'Rion', 'Rion_s', 'Cion_s', 'Cchem', 'Yeon_p', 'neon_p', 'R0'};
        lowerBounds    = [     0,        0,        0,       0,        0,      0.0,    0];
        upperBounds    = [   Inf,      Inf,      Inf,     Inf,      Inf,      1.0,  Inf];
    case '7c'
        pNames         = {'Rion', 'Rion_s', 'Yion_s', 'nion_s', 'Cchem', 'Ceon_p', 'R0'};
        lowerBounds    = [     0,        0,        0,        0,       0,      0.0,    0];
        upperBounds    = [   Inf,      Inf,      Inf,        1,     Inf,      Inf,  Inf];
    case 'maier2006'
        pNames         = {'R_{lyte}', 'R_{ion,int}', 'Q_{int}', 'n_{ion,int}', 'R_{ion,surf}', 'Q_{chem}', 'n_{chem}'};
        lowerBounds    = [         1,             1,         0,           0.0,              1,        0.0,        0.0];
        upperBounds    = [       Inf,           Inf,       Inf,           1.0,            Inf,        Inf,        1.0];
end

% Scale everything so it falls in the range 1:10
parameters = startingParams;
fitParamMatrix = zeros(length(runs), length(startingParams));       
scalingFactors = 10.^floor(log10(parameters));

%% Open summary filehandle for writing
fitsSummaryFilename = [sampleDir slash fitFolderName slash fitFolderName '.txt']
if(writeToFile)
    sumFH = fopen(fitsSummaryFilename, 'w');
    
    fprintf(sumFH, '"Filename"\t"Chi-Sqr"\t"Sum-Sqr"\t');
    for i = 1:length(pNames)
        head = pNames{i};
        fprintf(sumFH, '%s(+)\t%s(Error)\t%s(Error%%)\t', head, head, head);
    end
    fprintf(sumFH, '\r\n');
end
%% Do the fitting!
bar = waitbar(0, 'Initializing...');

for i = 1:length(runs)
    if runs(i).customFit == 1
        % Get frequency and impedance data in appropraite range
        simName = [sampleDir slash fitFolderName slash runs(i).filename '.sim'];
        freqs = runs(i).Z.freq;
        fRange = getRangeIndices(freqs, freqRange);
        %closestValues = min(abs(runs(i).Z.im(fRange))); %%%
        %indexOfClosest = find(closestValues == abs(runs(i).Z.im));
        %fRange = (max(1,indexOfClosest-1)):1:min(indexOfClosest+1,length(freqs)); %%%

        
        % Calculate measured and initial impedance profiles
        % But only include where Zim > 0
        
        posIndices = runs(i).Z.im(fRange) < 0;
        newFRange = [];
        for jj = 1:length(posIndices)
            if(posIndices(jj))
                newFRange = [newFRange; fRange(jj)];
            end
        end
        freqs = freqs(newFRange);
        
        Zmeas = [runs(i).Z.re(newFRange), runs(i).Z.im(newFRange)];
        Zinitial = Circuits(parameters, freqs, circuitVersion);
        
        % Scale parameters and bounds by appropriate factors
        scaledParameters = parameters ./ scalingFactors;
        scaledLowerBounds = lowerBounds ./ scalingFactors;
        scaledUpperBounds = upperBounds ./ scalingFactors;
        
        waitbar(i/length(runs), bar, sprintf('Fitting spectrum %d of %d...', i, length(runs)));
        
        % Weighted arc Fitting 'normalizes' the parameters by dividing by the scaling factor
        [p_out, resnorm, residuals, eflag, output, lambda, jacob] = ...
            lsqnonlin(@(ps) WeightedArcFitting(ps, scalingFactors, freqs, circuitVersion, Zmeas, weighting), ...
            scaledParameters, scaledLowerBounds, scaledUpperBounds, opts);
        fitParams = p_out .* scalingFactors;
        covar = (jacob' * jacob)^(-1);
        ci = nlparci(p_out, residuals, 'jacobian', jacob, 'alpha',  0.33);
        cLow = ci(:, 1)' .* scalingFactors;
        cHigh = ci(:, 2)' .* scalingFactors;
        %fittedFreqs = logspace(0, 8)';
        
        fittedFreqs = freqs;
        Zfitted = Circuits(fitParams, fittedFreqs, circuitVersion);
        fitParamMatrix(i, :) = fitParams;
        
        %% Plot the results
        if(makePlots)
            nyFH = figure();
            set(nyFH, 'color', 'white');
            axis square;
            set(gca, 'DataAspectRatio', [1 1 1]);
            
            hold on;
            plot(Zmeas(:, 1),    -Zmeas(:, 2), 'co');
            plot(Zinitial(:, 1), -Zinitial(:, 2), 'r--');
            plot(Zfitted(:, 1),  -Zfitted(:, 2), 'b-');
            hold off;
            legend('Measured', 'Guessed', 'Fitted', 'Location', 'SouthEast')
            title(num2str(runs(i).seqNum));
        end
        
        %% Show the results
        if(verbose)
            disp(['Run #' num2str(runs(i).seqNum)]);
            for k = 1:length(pNames)
                str = [pNames{k}, ': ', num2str(fitParams(k)), ', [', ...
                    num2str(cLow(k)), ', ' num2str(cHigh(k)) ']'];
                disp(str);
            end
            disp(' ');
        end
        
        %% Write to file
        if(writeToFile)
            fprintf(sumFH, '"%s"\t%E\t%E\t', ...
                [sampleDir slash fitFolderName slash runs(i).filename], 0, resnorm);
            for k = 1:length(fitParams)
                errorAbs = (cHigh(k) - cLow(k)) / 2;
                errorPercent = errorAbs / fitParams(k);
                fprintf(sumFH, '%E\t%E\t%E\t', fitParams(k), errorAbs, errorPercent);
            end
            fprintf(sumFH, newline);
            
            simFH = fopen(simName, 'w');
            for k = 1:10
                fprintf(simFH, ['Hey you! Keep on keepin'' on!' newline]);
            end
            fprintf(simFH, ['Freq(Hz)\tAmpl\tBias\tTime(Sec)\tZ''(a)\tZ''''(b)\tGD\tErr\tRange' newline]);
            for k = 1:length(fittedFreqs)
                fprintf(simFH, ['%E, 0.0E+00, 0.0E+00, 0.0E+00, %E, %E, 0.0E+00, 0, 0' newline], fittedFreqs(k), Zfitted(k, 1), Zfitted(k, 2));
            end
            fclose(simFH);
        end
    end
end
if(writeToFile)
    fclose(sumFH);
end
close(bar);

nonzeroEntries = find(fitParamMatrix(:, 1));
medianFitParams = median(fitParamMatrix(nonzeroEntries, :));

end

