% This script applies the given rotation matrix around the given point to the rows of the data matrix

function rotated_matrix = Rotation_Matrix_Application(rotation_matrix, rotation_point, data_matrix)

    %% Rotation %%
        % Translation w.r.t. the rotational point
        data_matrix     = data_matrix - rotation_point;
    
        % Matrix multiplication
        rotated_matrix  = (rotation_matrix * data_matrix')';
        
        % Translated back
        rotated_matrix  = rotated_matrix + rotation_point;

end