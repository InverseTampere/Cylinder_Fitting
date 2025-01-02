% The projections of vectors onto a plane are the vectors themselves minus the projections of the vector onto the plane's normal vector
% Note that the first output is not necessarily of unit length

function [proj_vector_matrix, proj_unit_vector_matrix, delta_list] = Vector_to_Plane_Projection(vector_matrix, plane_normal_vector)
    
    %% Projection onto plane %%
        % The projections onto the normal vector
        vector_plane_normals    = Vector_to_Vector_Projection(vector_matrix, plane_normal_vector);

        % The resulting projection onto the plane
        proj_vector_matrix      = vector_matrix - vector_plane_normals;
        delta_list              = sqrt(sum(proj_vector_matrix.^2, 2));
        proj_unit_vector_matrix = proj_vector_matrix ./ delta_list;
        
end