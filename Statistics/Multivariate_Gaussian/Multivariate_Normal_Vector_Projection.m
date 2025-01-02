% The projection of an n-D MVN onto the given vectors running through the given point is performed using Pope's method [O11]
% Note that the minimum and maximum extrema are given in the third dimension of the matrix

function [projected_centre_matrix, projected_extrema_matrix, projection_extent_list] = Multivariate_Normal_Vector_Projection(vector_matrix, point, mu, covariance_matrix)

    %% Cholesky factorisation %%
        % The lower-diagonal Cholesky factorisation of the ellipsoidal matrix
        Ellipsoid_matrix    = inv(covariance_matrix);
        L                   = chol(Ellipsoid_matrix, 'lower')';   
        
    %% Extent of the ellipsoid onto the vectors %%
        % The vectors are normalised
        vector_matrix = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));
        
        % Projection of the ellipsoid onto the vectors
        proj_ellipsoid_vectors  = (L \ vector_matrix')';
        
        % The extent of the projected ellipsoid is twice the norms of these vectors
        projection_extent_list  = 2 * sqrt(sum(proj_ellipsoid_vectors.^2, 2));
        
    %% Ellipsoid centre projection %%
        % Projection of the ellipsoid centre onto the vectors
        Centre_Projection_fun   = @(vector) Point_to_Vector_Projection(mu, vector, point);
        
        [num_vectors, num_dim]  = size(vector_matrix);
        vector_cell             = mat2cell(vector_matrix, ones(1, num_vectors), num_dim);
        projected_centre_matrix = cellfun(Centre_Projection_fun, vector_cell, 'UniformOutput', false);
        projected_centre_matrix = cell2mat(projected_centre_matrix);
        
    %% Extrema of projected ellipsoid %%
        % The resulting extrema on the vectors
        projected_extrema_matrix = projected_centre_matrix + vector_matrix .* cat(3, -projection_extent_list, projection_extent_list)/2; 
        
end