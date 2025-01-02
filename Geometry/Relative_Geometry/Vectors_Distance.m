% This script computes the distance between two (sets of) vectors or between a vector and a point
% In the latter case, the vector of said point can be left as a zero-vector

% It does this by creating the vector that is orthogonal to both, and gives the intersections of this spanning vector to both given vectors

function [distance_list, inter_a_matrix, lambda_a_list, inter_b_matrix, delta_b_list] = Vectors_Distance(vector_a_matrix, a_start_matrix, vector_b_matrix, b_start_matrix)
    
    %% Vectors are normalised %%
        % Division by the vector norms
        vector_a_matrix = vector_a_matrix ./ sqrt(sum(vector_a_matrix.^2, 2));
        vector_b_matrix = vector_b_matrix ./ sqrt(sum(vector_b_matrix.^2, 2));
        
        % Zero-vectors create NaN values. These are replaced by zeros again
        vector_a_matrix(isnan(vector_a_matrix)) = 0;
        vector_b_matrix(isnan(vector_b_matrix)) = 0;

    %% Intersection points on vectors a and b %%
        % The distances along vectors a and b at which the shorted distance between them occurs
        delta_b_list    = (dot(b_start_matrix - a_start_matrix, vector_a_matrix, 2) .* dot(vector_a_matrix, vector_b_matrix, 2) - dot(b_start_matrix - a_start_matrix, vector_b_matrix, 2)) ./ (1 - dot(vector_b_matrix, vector_a_matrix, 2).^2);      % Distance along vector b
        lambda_a_list   = dot(b_start_matrix - a_start_matrix, vector_a_matrix, 2) + delta_b_list .* dot(vector_a_matrix, vector_b_matrix, 2);         % Distance along vector a
                    
        % The points where the shortest vector intersects vectors a and b
        inter_a_matrix  = a_start_matrix + lambda_a_list .* vector_a_matrix;
        inter_b_matrix  = b_start_matrix + delta_b_list .* vector_b_matrix;
    
    %% Shortest distance between vectors a and b %%
        % Norm of the vector between the intersects
        shortest_vector_matrix  = inter_a_matrix - inter_b_matrix;
        distance_list           = sqrt(sum(shortest_vector_matrix.^2, 2));       

end