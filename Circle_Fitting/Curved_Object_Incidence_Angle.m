% This script computes the incidence angle between the given vectors and a circle at its intersection point
% The output incidence angles are 0 if the vector is orthogonal to the surface

function incidence_angle_list = Curved_Object_Incidence_Angle(x_i_list, y_i_list, vector_x_list, vector_y_list, circle_x, circle_y)
    %% Column vectors are ensured %%
        number_vectors  = length(vector_x_list);

        x_i_list        = reshape(x_i_list, [number_vectors, 1]);
        y_i_list        = reshape(y_i_list, [number_vectors, 1]);
        vector_x_list   = reshape(vector_x_list, [number_vectors, 1]);
        vector_y_list   = reshape(vector_y_list, [number_vectors, 1]);

    %% Normalised vectors %%
        % The normalised vectors between the circle centre and intersection point
        circle_intersection_vectors = [x_i_list - circle_x, y_i_list - circle_y];
        circle_intersection_vectors = circle_intersection_vectors ./ sqrt(sum(circle_intersection_vectors.^2, 2));

        % The given vectors are also normalised
        vectors_matrix = [vector_x_list, vector_y_list];
        vectors_matrix = vectors_matrix ./ sqrt(sum(vectors_matrix.^2, 2));
    
    %% The incidence angles %%
        % The smallest angle between the two vectors is the incidence angle
        incidence_angle_list = acos(abs(circle_intersection_vectors * vectors_matrix'));
        incidence_angle_list = diag(incidence_angle_list);

        incidence_angle_list = real(incidence_angle_list);      % Sometimes a negligible imaginary component can be present
end