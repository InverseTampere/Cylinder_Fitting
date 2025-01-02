% The expected Mahalanobis distance over the given circle is computed numerically by sampling the distribution
% To increase computational speed, the circle may also be sampled if Circle_Sampling is true

function expected_Mahal_distance = Expected_Mahalanobis_Distance_Circle_Numerical(circle_centre, circle_radius, distr_mu, distr_sigmae, distr_axes, distance_moment, number_samples, Confidence_interval, Circle_Sampling, Diagnostics)
    
    %% Manual inputs %%
        % The number of samples placed on the circle can be less than the given number of samples
        number_samples_circle = 1e2;

    %% Coordinate frame transformation %%
        % The Mahalanobis distance equals the Euclidean distance when the distribution is transformed into a standard-normal
        num_dim             = length(distr_mu);
        squeeze_matrix      = 1./distr_sigmae .* eye(num_dim);
    
        % The circle (now ellipse)'s properties
        circle_centre_t     = circle_centre - distr_mu;
        circle_centre_r     = (distr_axes * circle_centre_t')';
        ellipse_centre      = (squeeze_matrix * circle_centre_r')';
        
        ellipse_axes        = eye(num_dim);
        ellipse_radii       = circle_radius ./ distr_sigmae;
        
    %% Standard normal sampling %%
        % Confidence interval
        m_STD = sqrt(chi2inv(Confidence_interval/100, num_dim));
        
        % Sampling the standard-normal distribution placed at the origin
        origin                          = zeros(1, num_dim);
        [sample_coord_matrix, ~, dA]    = Equal_Area_Circular_Sampler(m_STD, origin, [], number_samples);
        
        % The probability density of each sample
        prob_dens_list = 1/(2*pi) * exp(-1/2*sum(sample_coord_matrix.^2, 2));
                
    %% Expected Mahalanobis distance %%
        % Euclidean distance from each sample to the ellipse (previously circle), which equals the Mahalanobis distance in the original frame
        if Circle_Sampling == true
            % The distance is estimated numerically to avoid solving a fourth order polynomial
            [~, Mahal_dist_list] = Point_to_Ellipse_Projection_Numerical(ellipse_centre, ellipse_radii, ellipse_axes, sample_coord_matrix(:, 1:2), number_samples_circle, Diagnostics);
        else
            % The distance is determined analytically
            Projection_fun          = @(x, y) Point_to_Ellipse_Projection(ellipse_centre, ellipse_radii, ellipse_axes, [x, y], Diagnostics, Diagnostics);
            [~, Mahal_dist_cell]    = arrayfun(Projection_fun, sample_coord_matrix(:, 1), sample_coord_matrix(:, 2), 'UniformOutput', false);
            Mahal_dist_list         = vertcat(Mahal_dist_cell{:});
        end
        
        % Expected Mahalanobis distance
        expected_Mahal_distance = sum(prob_dens_list .* Mahal_dist_list.^distance_moment) * 100 / Confidence_interval * dA;         
                
end