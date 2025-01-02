% This script computes the normal vectors to the given triangles

function normal_vector_matrix = Triangle_Normal_Vector(Triangle_matrix, coordinate_matrix)

    %% Surface normal vectors %%
        % Each triangle consists of three vectors, any two of which can be used to compute its normal vector
        
        % The triangle vertices
        coord_a_matrix  = coordinate_matrix(Triangle_matrix(:, 1), :);
        coord_b_matrix  = coordinate_matrix(Triangle_matrix(:, 2), :);
        coord_c_matrix  = coordinate_matrix(Triangle_matrix(:, 3), :);
        
        % Two triangle vectors per triangle
        vector_ab_matrix = coord_b_matrix - coord_a_matrix;      % Vectors from a to b
        vector_ac_matrix = coord_c_matrix - coord_a_matrix;      % Vectors from a to c

        % The normal vectors are computed through the cross-product
        normal_vector_matrix = cross(vector_ab_matrix, vector_ac_matrix, 2);
        normal_vector_matrix = normal_vector_matrix ./ sqrt(sum(normal_vector_matrix.^2, 2));       % Normalised
                        
end