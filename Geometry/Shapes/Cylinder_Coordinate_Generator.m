% This script creates coordinates in (x, y, z) of the given cylinder, discretised in circles spanning its length

function [cylinder_x_list, cylinder_y_list, cylinder_z_list] = Cylinder_Coordinate_Generator(cyl_centre, cyl_radius, cyl_length, cyl_direction, number_circle_steps, number_length_steps)

    %% Cylinder coordinates %%
        % The coordinates are generated for each circle along the length of the cylinder
        theta_list      = linspace(0, 2*pi, number_circle_steps + 1)';
        theta_list(number_circle_steps + 1) = [];       % To avoid duplicate points at 0 and 2pi, the final point is removed
        
        length_steps    = linspace(-cyl_length/2, cyl_length/2, number_length_steps);
        
        circle_x_cell   = cell(1, number_length_steps);
        circle_y_cell   = cell(1, number_length_steps);
        circle_z_cell   = cell(1, number_length_steps);
        
        for l = 1 : number_length_steps
            length = length_steps(l);
            
            % The circle in the (x, y) plane
            circle_coordinates = [cyl_radius * cos(theta_list), cyl_radius * sin(theta_list), length * ones(number_circle_steps, 1)];
            
            % And rotated to match the cylinder direction
            vector_start            = zeros(1, 3);     % Note that the circle coordinates are defined w.r.t. the origin
            coordinate_matrix_rot   = Rotation_3D(circle_coordinates, cyl_direction, vector_start);
            
            % The centre of the cylinder is added to the final coordinates
            circle_x_cell{l} = coordinate_matrix_rot(:, 1) + cyl_centre(1);
            circle_y_cell{l} = coordinate_matrix_rot(:, 2) + cyl_centre(2);
            circle_z_cell{l} = coordinate_matrix_rot(:, 3) + cyl_centre(3);
        end
    
        % The coordinates are merged
        cylinder_x_list = vertcat(circle_x_cell{:});
        cylinder_y_list = vertcat(circle_y_cell{:});
        cylinder_z_list = vertcat(circle_z_cell{:});
        
end