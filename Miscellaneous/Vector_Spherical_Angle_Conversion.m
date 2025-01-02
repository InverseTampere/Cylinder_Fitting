% This script converts between a given 3D vector and the azimuth and elevation angles
% The given vector does not need to be unit length, but that a unit vector is always produced from spherical angles
% When converting from vector to spherical angles, only the vector needs to be given as input and vice versa

% Note: The non-signed (ns) outputs are sign-independent of the true direction (i.e. they may point in exactly the opposite direction)
%       For the non-signed (ns) outputs the elevation angle is in range [0, pi/2] and azimuth angle in range [-pi, pi]. The output vector therefore always has positive z
%       For the signed (s) outputs the elevation angle is in range [-pi/2, pi/2] and azimuth angle in range [-pi, pi]. The output vector may be of any sign
    
function [vector_matrix_ns, azim_angle_list_ns, elev_angle_list_ns, vector_matrix_s, azim_angle_list_s, elev_angle_list_s] = Vector_Spherical_Angle_Conversion(vector_matrix, azim_angle_list, elev_angle_list)

    %% Vector to spherical angle conversion %%
        if ~isempty(vector_matrix)
            % Converted to unit length
            vector_matrix               = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));          

            % Non-signed outputs
            sign_list                   = sign(vector_matrix(:, 3));                                % Third dimension is ensured to be positive
            sign_list(sign_list == 0)   = 1;                                                        % In case the z-component is zero, the vectors are unchanged
            vector_matrix_ns            = sign_list .* vector_matrix;                               

            azim_angle_list_ns  = atan2(vector_matrix_ns(:, 2), vector_matrix_ns(:, 1));            % [-pi, pi] 
            elev_angle_list_ns  = asin(vector_matrix_ns(:, 3));                                     % [0, pi/2]

            % Signed outputs
            azim_angle_list_s   = atan2(vector_matrix(:, 2), vector_matrix(:, 1));                  % [-pi, pi]
            elev_angle_list_s   = asin(vector_matrix(:, 3));                                        % [-pi/2, pi/2]

            vector_matrix_s     = vector_matrix;
        end

    %% Spherical angle to vector conversion %%
        if ~isempty(azim_angle_list) && ~isempty(elev_angle_list)
            % Signed outputs
            vector_matrix_s     = [cos(elev_angle_list) .* cos(azim_angle_list), ...
                                   cos(elev_angle_list) .* sin(azim_angle_list), ...
                                   sin(elev_angle_list)];

            elev_angle_list_s   = elev_angle_list;
            azim_angle_list_s   = azim_angle_list;

            % Non-signed outputs
            sign_list                   = sign(elev_angle_list);                                    % The elevation angle is ensured to be positive
            sign_list(sign_list == 0)   = 1;                                                        % In case the elevation angle is zero, the angles are unchanged
            elev_angle_list_ns          = sign_list .* elev_angle_list;

            vector_matrix_ns            = sign_list .* vector_matrix_s;
            azim_angle_list_ns          = atan2(vector_matrix_ns(:, 2), vector_matrix_ns(:, 1));
        end        
end
