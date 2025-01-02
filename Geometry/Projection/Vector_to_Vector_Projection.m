% This script gives the projection of a set of projected vectors onto another projection vector

function projected_vector_matrix = Vector_to_Vector_Projection(vector_matrix, projection_vector)

    %% Vector projection %%
        % For efficiency, the projection vector is repeated
        num_vectors             = size(vector_matrix, 1);
        projection_vector_rep   = repmat(projection_vector, [num_vectors, 1]);
                
        % The projected vectors
        projected_vector_matrix = dot(vector_matrix, projection_vector_rep, 2)/norm(projection_vector)^2 .* projection_vector;
        
end