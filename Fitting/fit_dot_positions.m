function [ residuals ] = fit_dot_positions( fit_parameters, indices, positions )
%FIT_DOT_POSITIONS Wrapper for LSQNONLIN to fit the measured dot positions
%to  a global model
%   Inputs:
%       fit_parameters - array of fit parameters to pass to calculate
%       transformation
%       indices - matrix of dot indices for which postiions were measured
%       positions - the positions for each measured dot
%   Outputs:
%       residuals - the difference between fitted and actual dot position


% Compute the transformation of all dots based on the given parameters
dot_position_matrix = calculate_dot_positions(fit_parameters);

% Compute the residual distance between the transformed dots and the measured dots
[nDots, ~] = size(indices);
residuals = zeros(3, nDots);
for dot = 1:nDots
    row = indices(dot, 1);
    column = indices(dot, 2);
    residuals(:, dot) = squeeze(dot_position_matrix(row, column, :)) - positions(:, dot);
end

end

