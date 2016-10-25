function [resVec] = ...
    WeightedArcFitting(parameters, scalingFactors, freqs, ...
                        circuitVersion, Zmeas, weight)
%WEIGHTEDARCFITTING returns residuals of calculated impedance
%   To improve the fitting ability of lsqnonlin, it is helpful to break up
%   the data into real and imaginary residuals for each point. This
%   function determines these residuals between measured data (Zmeas) and
%   the calculated data (from Circuits). The user should create an
%   anonymous function to capture the output of this function given the
%   options selected by the user.

parameters = parameters .* scalingFactors;

Zcalc = Circuits(parameters, freqs, circuitVersion);
resVec = zeros(length(freqs), 2);

switch weight
    case 'modulus'  % modulus weighting
        resVec = (Zmeas - Zcalc) ./ repmat(sqrt(sum(Zmeas'.^2)), 2, 1)';
        %for k = 1:length(freqs)
        %    denom = sqrt(sum(Zmeas(k, :).^2));
        %    resVec(2*k - 1) = (Zmeas(k,1) - Zcalc(k,1)) / denom;
        %    resVec(2*k)     = (Zmeas(k,2) - Zcalc(k,2)) / denom;
        %end
    case 'unit'  % unit weighting
        for k = 1:length(freqs)
            resVec(2*k - 1) = (Zmeas(k,1) - Zcalc(k,1));
            resVec(2*k)     = (Zmeas(k,2) - Zcalc(k,2));
        end
    case 'proportional' % proportional weighting
        for k = 1:length(freqs)
            resVec(2*k - 1) = (Zmeas(k,1) - Zcalc(k,1)) ./ Zcalc(k,1);
            resVec(2*k)     = (Zmeas(k,2) - Zcalc(k,2)) ./ Zcalc(k,2);
        end
    otherwise
        errstr = ['Unrecognized weight string "' weight '.'];
        error(errstr);
end

end