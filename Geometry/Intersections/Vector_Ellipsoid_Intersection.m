% The intersection between a vector starting from a point and an ellipsoid can be computed in a straightforward manner, once the ellipsoid has been transformed into a unit sphere
% The first and second intersections are given in the 3rd dimension of the intersection matrix respectively
% If the intersection does not exist, NaN values are given instead

function [intersection_matrix, delta_matrix] = Vector_Ellipsoid_Intersection(ellipsoid_centre, ellipsoid_axes, ellipsoid_radii, vector, point_matrix)

    %% Transformation matrices %%    
        % Rotation s.t. the vector forms the z-axis
        vector = vector / norm(vector);
        z_v_vector = vector;
        
        projected_ellipsoid_axes = Vector_to_Plane_Projection(ellipsoid_axes, vector);
        x_v_vector = projected_ellipsoid_axes(1, :);
        x_v_vector = x_v_vector / norm(x_v_vector);

        y_v_vector = cross(z_v_vector, x_v_vector);

        vector_basis = [x_v_vector; y_v_vector; z_v_vector];
        
        % Matrix to squeeze to a unit sphere
        [num_samples, num_dim]  = size(point_matrix);
        squeeze_matrix          = 1./ellipsoid_radii .* eye(num_dim);
        
        vector_s                = (squeeze_matrix * ellipsoid_axes * vector')';
                
        % The full transformation matrix
        transformation_matrix = vector_basis * squeeze_matrix * ellipsoid_axes;
                
    %% Intersections %%
        % For the integrals to be nonzero, the samples must be located within the unit sphere
        % The distance (delta) from the point along the vector to the intersections with the unit sphere is calculated in the vector aligned frame
        
        % The points and vector are moved to the vector-aligned frame
        vector_v        = (transformation_matrix * vector')';
        point_matrix_t  = point_matrix - ellipsoid_centre;
        point_matrix_v  = (transformation_matrix * point_matrix_t')';
                    
        % The intersections with the unit sphere
        unit_sphere_radius  = 1;
        origin              = zeros(1, num_dim);
        vector_matrix_v     = repmat(vector_v, [num_samples, 1]);

        [~, lambda_matrix]  = Vector_Sphere_Intersection(point_matrix_v, vector_matrix_v, origin, unit_sphere_radius);
        
        % The distances are compensated for squeezing
        delta_matrix        = lambda_matrix / norm(vector_s);
    
        % The resulting intersections
        lambda_matrix_res   = reshape(delta_matrix, [num_samples, 1, 2]);
        intersection_matrix = point_matrix + vector .* lambda_matrix_res;
        
end