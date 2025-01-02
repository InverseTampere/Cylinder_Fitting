% Vector bases are computed for the given cylinder and scanner locations
% The scanner vector bases are such that the scanner lies on the y-z plane
% The spectral angles of the cylinder axis are given in this coordinate frame
% The projection vector bases are the full rotation matrices to the cross-sectional coordinate frame

function [scanner_vector_base_cell, azim_angle_list, elev_angle_list, projection_vector_base_cell] = Vector_Basis_Cylinder_Cross_Section_Projection(cylinder_centre, cylinder_direction, Scanner_loc_cell)

    %% Manual inputs %%
        numerical_margin    = 1e-6;     % [-] To check for gimbal lock

    %% Vector bases for each scanner %%
        % Normalising the input vector
        cylinder_direction = cylinder_direction / norm(cylinder_direction);

        % Properties for each scanner
        number_scanners = length(Scanner_loc_cell);

        scanner_vector_base_cell    = cell(1, number_scanners);
        projection_vector_base_cell = cell(1, number_scanners);
        azim_angle_list             = zeros(1, number_scanners);
        elev_angle_list             = zeros(1, number_scanners);

        for s = 1 : number_scanners
            % Relative scanner location w.r.t. the cylinder centre
            scanner_loc     = Scanner_loc_cell{s} - cylinder_centre;
            [s_i, s_j, s_k] = Column_Deal(scanner_loc);

            % Rotation matrix s.t. the scanner lies on the y-z plane    
            if s_j == 0 && s_i ~= 0
                scanner_vector_basis = sign(s_i) * [0, 1, 0; 1, 0, 0; 0, 0, 1];             % First and second dimensions are permuted and potentially all axes have their sign changed
            elseif s_i == 0 && s_j ~= 0
                scanner_vector_basis = [sign(s_j), 0, 0; 0, sign(s_j), 0; 0, 0, 1];         % No rotation, but potentially the first and second dimensions have their signs changed
            elseif s_i == 0 && s_j == 0
                scanner_vector_basis = sign(s_k) * [1, 0, 0; 0, 0, 1; 0, 1, 0];             % Second and third dimensions are permuted and potentially all axes have their sign changed
            else
                a_i     = s_j / sqrt(s_i^2 + s_j^2);
                a_j     = -sqrt(1 - a_i^2);

                scanner_vector_basis = [a_i, a_j, 0; -a_j, a_i, 0; 0 ,0, 1];
            end

            scanner_vector_base_cell{s} = scanner_vector_basis;

            % Rotating the cylinder axis and computing the spherical angles
            cyl_dir_0 = (scanner_vector_basis * cylinder_direction')';
    
            [~, ~, ~, ~, gamma, epsilon]    = Vector_Spherical_Angle_Conversion(cyl_dir_0, [], []);
            azim_angle_list(s)              = gamma;
            elev_angle_list(s)              = epsilon;
            
            % Vector basis s.t. the scanner still lies in the y-z plane, with the z-axis being parallel to the cylinder axis
            cyl_dir_a = cyl_dir_0(1);

            if abs(cyl_dir_a) > 1 - numerical_margin            % If the cylinder direction is parallel to the a-axis, the equations below do not work
                x_vec   = [0, 0, -sign(cyl_dir_a)];
                y_vec   = [0, 1, 0];
            else
                nu      = 1/sqrt(sin(epsilon)^2 + cos(epsilon)^2*sin(gamma)^2);
                x_vec   = nu * [sin(epsilon)^2 + cos(epsilon)^2*sin(gamma)^2, -cos(epsilon)^2*cos(gamma)*sin(gamma), -cos(epsilon)*sin(epsilon)*cos(gamma)];
                y_vec   = nu * [0, sin(epsilon), -cos(epsilon)*sin(gamma)];        
            end

            % Cylinder vector basis
            cylinder_vector_basis = [x_vec; y_vec; cyl_dir_0];    

            % Full projection vector_basis
            projection_vector_basis         = cylinder_vector_basis * scanner_vector_basis;
            projection_vector_base_cell{s}  = projection_vector_basis;
        end
end