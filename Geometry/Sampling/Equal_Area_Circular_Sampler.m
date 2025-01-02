% This script creates samples that are equally spaced and covering a circular area in 2D or 3D (requiring a normal vector)
% The total number of samples may differ slightly from the initially given value

function [sample_coord_matrix, number_samples, delta_area] = Equal_Area_Circular_Sampler(circle_radius, circle_centroid, normal_vector, number_samples)
        
    %% Circle on x-y plane %%
        % The number of samples in each circle is equivalent to a square, with the following samples per side of each quarter
        n_quarter = floor(2/pi * sqrt(number_samples));

        % The central row and column
        var_list    = circle_radius/(n_quarter + 1) * linspace(1, n_quarter + 1, n_quarter);
        
        x_list_c    = [-var_list, 0, var_list];         % Note how the origin is added as well
        y_list_c    = [-var_list, 0, var_list];
        
        % The full grid
        [x_matrix, y_matrix] = meshgrid(x_list_c, y_list_c);
        
        % Of all the samples, only the ones with norm within the radius are kept
        norm_matrix     = sqrt(x_matrix.^2 + y_matrix.^2);
        
        x_circle_list   = x_matrix(norm_matrix <= circle_radius);
        y_circle_list   = y_matrix(norm_matrix <= circle_radius);
        
        number_samples  = length(x_circle_list);
        
        % The area per sample is then
        delta_area = pi * circle_radius^2 / number_samples;
        
    %% Rotation and translation %%
        % Translated according to the centroid
        sample_coord_matrix = [x_circle_list, y_circle_list];
        sample_coord_matrix = sample_coord_matrix + circle_centroid;

        % The samples are rotated s.t. the normal vector forms the z-axis
        if ~isempty(normal_vector)
            % A third dimension of zeros is added
            sample_coord_matrix = [sample_coord_matrix, zeros(number_samples, 1)];

            % They are rotated        
            [sample_coord_matrix, ~] = Rotation_3D(sample_coord_matrix, normal_vector, circle_centroid);
        end
end