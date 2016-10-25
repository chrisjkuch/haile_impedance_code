function [ positions ] = calculate_dot_positions( transformation_parameters )
%CALCULATE_DOT_POSITIONS Calculate the dot positions
%   Inputs:
%       alpha - coefficient of thermal expansion
%       theta_zz - rotation about the z-axis
%       theta_yy - rotation about the x-axis
%       theta_xx - rotation about the y-axis
%       s_x, _y, _z - substrate position
%   Outputs:
%       positions - matrix of dot positions

alpha_cte = transformation_parameters(1);
theta_zz = transformation_parameters(2);
theta_yy = transformation_parameters(3);
theta_xx = transformation_parameters(4);
substrate_x = transformation_parameters(5);
substrate_y = transformation_parameters(6);
substrate_z = transformation_parameters(7);
offset = [substrate_x; substrate_y; substrate_z];
% Initialize all dot positions
x_pos = -[1 1.85:0.350:8.25 9];
y_pos = -[0.8 1.3 1.7 2 2.3 2.5:0.2:3.1 3.25:0.15:4.3];

nRows = length(y_pos);
nColumns = length(x_pos);
dot_positions = zeros(nRows, nColumns, 3);

% Generate rotation matrices
R_Z = [[cos(theta_zz), -sin(theta_zz), 0]; [sin(theta_zz), cos(theta_zz), 0]; [0 0 1]];
R_Y = [[cos(theta_yy), 0, sin(theta_yy)]; [0, 1, 0]; [-sin(theta_yy), 0, cos(theta_yy)]];
R_X = [[1, 0, 0]; [0, cos(theta_xx), -sin(theta_xx)]; [0, sin(theta_xx), cos(theta_xx)]];

% Transform all of the points
for row = 1:nRows
    for col = 1:nColumns
        position_vector = alpha_cte * [x_pos(col); y_pos(row); 0];
        position_vector = R_Z * R_Y * R_X * position_vector + offset;
        dot_positions(row, col, :) = position_vector;
    end
end

positions = dot_positions;

end