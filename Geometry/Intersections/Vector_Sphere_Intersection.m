% Between a vector and a sphere there can be zero, one or two intersections which are returned by this script
% They are given in the 3rd dimension of the intersection points matrix
% Non-intersects are given by NaN
% Lambda is the distance along the vector to the intersections

function [intersect_points_matrix, first_intersect_point_matrix, lambda_matrix, intersection_bool] = Vector_Sphere_Intersection(vector_start_matrix, vector_matrix, sphere_centre, sphere_radius)
    
    %% Intersection points %%
        % The vector is normalised
        vector_matrix   = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));
    
        % The distance along the vector to the intersections
        lambda_matrix    = -dot(vector_matrix, vector_start_matrix - sphere_centre, 2) + [-1, 1] .* sqrt(dot(vector_matrix, vector_start_matrix - sphere_centre, 2).^2 - sum((vector_start_matrix - sphere_centre).^2, 2) + sphere_radius^2);
        
        % Only real delta's can lead to a valid intersection
        valid_deltas                    = imag(lambda_matrix) == 0;
        lambda_matrix(~valid_deltas)    = NaN;

        intersection_bool               = sum(valid_deltas, 2) > 0;
    
        % The resulting valid intersections
        num_samples             = size(vector_start_matrix, 1);
        delta_matrix_res        = reshape(lambda_matrix, [num_samples, 1, 2]);
        intersect_points_matrix = vector_start_matrix + delta_matrix_res .* vector_matrix;

        first_intersect_point_matrix = squeeze(intersect_points_matrix(:, :, 1));
end