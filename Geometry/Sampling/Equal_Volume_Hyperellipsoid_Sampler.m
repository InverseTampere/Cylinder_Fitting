% This script creates samples that are equally spaced and covering an n-dimensional hyperellipsoidal volume
% Note that the final number of samples usually differs slightly from the original input

function [ellipsoid_sample_matrix, number_samples, delta_volume] = Equal_Volume_Hyperellipsoid_Sampler(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, number_samples)
        
    %% Sampling enclosing hypercuboid %%
        % Volume of the cuboid
        cuboid_vertex_lengths   = 2*ellipsoid_radii;
        V_cuboid                = prod(cuboid_vertex_lengths);

        % Volume of the ellipsoid
        num_dim     = length(ellipsoid_centre);
        V_ellipsoid = pi^(num_dim / 2) / gamma(num_dim/2 + 1) * prod(ellipsoid_radii);

        % Number of samples required by the cuboid
        number_samples_cuboid   = V_cuboid / V_ellipsoid * number_samples;

        % Sampling the hypercuboid placed at the origin and unrotated
        origin                          = zeros(1, num_dim);
        identity_matrix                 = eye(num_dim);
        [cuboid_sample_matrix, ~, ~, ~] = Equal_Volume_Hypercuboid_Sampler(origin, cuboid_vertex_lengths, identity_matrix, number_samples_cuboid, []);
    
    %% Rejecting samples outside hyperellipsoid %%
        % Samples within the ellipse have an ellipsoidal value <= 1
        ellipsoidal_eq_matrix   = sum(cuboid_sample_matrix.^2 ./ ellipsoid_radii.^2, 2);
        ellipsoidal_boolean     = ellipsoidal_eq_matrix <= 1;

        ellipsoid_sample_matrix = cuboid_sample_matrix(ellipsoidal_boolean, :);

        % Number of samples remaining and their respective volume
        number_samples          = sum(ellipsoidal_boolean);
        delta_volume            = V_ellipsoid / number_samples;

    %% Transformation to original coordinate frame %%
        % First rotated, where the rotational matrix is the inverse of the ellipsoidal axes
        ellipsoid_sample_matrix = (ellipsoid_axes \ ellipsoid_sample_matrix')';

        % Then translated
        ellipsoid_sample_matrix = ellipsoid_sample_matrix + ellipsoid_centre;
        
end