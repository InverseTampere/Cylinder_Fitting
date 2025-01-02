% Rotation of a covariance matrix by the given rotation matrix

function Sigma_rotated = Covariance_Matrix_Rotation(Sigma, rotation_matrix)
    % As the covariance matrix is not simply in [x, y, z] format, the rotation matrix is applied twice
    Sigma_rotated = rotation_matrix * Sigma * rotation_matrix';             % Note that the transpose and inverse of a rotation matix are equal
end