% Euler rotation around a point is applied to the given coordinate matrix according to the provided angles
% Note that the angles are expected to be given in radians

function [coordinate_matrix_rot, rotation_matrix] = Euler_Rotation_3D(coordinate_matrix, rotation_point, pitch_angle, yaw_angle, roll_angle)

    %% Rotation matrix %%
        % Rotation components
        roll_matrix =   [1,                 0,                  0;
                         0,                 cos(roll_angle),    sin(roll_angle);
                         0,                 -sin(roll_angle),   cos(roll_angle)];

        pitch_matrix =  [cos(pitch_angle),  0,                  -sin(pitch_angle);
                         0,                 1,                  0;
                         sin(pitch_angle),  0,                  cos(pitch_angle)];

        yaw_matrix =    [cos(yaw_angle),    sin(yaw_angle),     0;
                         -sin(yaw_angle),   cos(yaw_angle),     0;
                         0,                 0,                  1];

        % The combined rotation matrix
        rotation_matrix = roll_matrix*pitch_matrix*yaw_matrix;

    %% Rotation of the coordinate matrix %%
        % Translation w.r.t. the rotational point
        coordinate_matrix = coordinate_matrix - rotation_point;
    
        % The function is applied to a cell array for efficiency
        rotation_fun = @(x) (rotation_matrix * x')';

        num_points              = size(coordinate_matrix, 1);
        coordinate_cell         = mat2cell(coordinate_matrix, ones(1, num_points));
        coordinate_cell_rot     = cellfun(rotation_fun, coordinate_cell, 'UniformOutput', false);
        coordinate_matrix_rot   = cell2mat(coordinate_cell_rot);
        
        % Translated back
        coordinate_matrix_rot   = coordinate_matrix_rot + rotation_point;
    
end
    