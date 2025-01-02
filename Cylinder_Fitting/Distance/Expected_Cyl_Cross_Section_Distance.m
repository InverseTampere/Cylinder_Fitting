% The expected distance from the projected ellipsoids to the (circular) cross-section
% The distance function can be complicated. Sampling is thus used

function [expected_distance, expected_distance_list, norm_expected_distance, norm_expected_distance_list] = Expected_Cyl_Cross_Section_Distance(Point_Cloud_Distributions, Fitting_Parameters, cylinder_centre, cylinder_radius, cylinder_direction, alpha, Print, Plot, Diagnostics)
        
    %% Inputs %%
    
        %--% Point cloud distributions %--%
        distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
        distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
        distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;
        distribution_Sigma_cell     = Point_Cloud_Distributions.distribution_Sigma_cell;
        number_distributions        = Point_Cloud_Distributions.num_distributions;
        
        %--% Fitting parameters %--%
        distance_moment             = Fitting_Parameters.distance_moment;
        number_samples              = Fitting_Parameters.num_ED_samples;
        ED_Bias_Correction          = Fitting_Parameters.ED_Bias_Correction;
        
    %% Translation according to the cylinder centre
        Translation_fun         = @(distribution_mu) distribution_mu - cylinder_centre;
        distribution_mu_c_cell  = cellfun(Translation_fun, distribution_mu_cell, 'UniformOutput', false);
        
    %% Projection onto the plane %%
        % Vector basis given the cylinder axis
        num_dim = length(cylinder_centre);
        origin  = zeros(1, num_dim);
        [~, vector_basis, ~] = Vector_Based_Rotation(origin, cylinder_direction, origin);
    
        % Projection of the distributions
        Distribution_Projection_fun     = @(distribution_mu, distribution_axes, distribution_sigmae, distribution_Sigma) Multivariate_Normal_Projection(vector_basis, origin, {distribution_mu}, {distribution_axes}, {distribution_sigmae}, {distribution_Sigma}, 1, Diagnostics);
        Projected_Distributions_cell    = cellfun(Distribution_Projection_fun, distribution_mu_c_cell, distribution_axes_cell, distribution_sigmae_cell, distribution_Sigma_cell, 'UniformOutput', false);

    %% Sampling the projected distributions %%
        % The number of standard deviations of the 2D ellipses given the confidence level
        m_STD = sqrt(chi2inv(1 - alpha, 2));
        
        % The resulting samples. Note that the amount may change slightly
        Sampling_fun        = @(Projected_Distributions) Equal_Area_Elliptical_Sampler(Projected_Distributions.Plane.Original.mu{1}, m_STD * Projected_Distributions.Plane.Original.sigmae{1}, Projected_Distributions.Plane.Original.distr_axes{1}, number_samples);       % Note conversion from cell array to regular arrays
        [samples_cell, ~]   = cellfun(Sampling_fun, Projected_Distributions_cell, 'UniformOutput', false);
        
    %% Integrated density values of the samples %%
        % The distributions are integrated along the cylinder axis
        Integration_fun             = @(samples_matrix, Projected_Distributions, distribution_Sigma, distribution_mu, distribution_axes, distribution_sigmae) Multivariate_Normal_Integration(samples_matrix, cylinder_direction, Projected_Distributions, distribution_Sigma, distribution_mu, distribution_axes, distribution_sigmae, Diagnostics);
        Statistical_Outputs_cell    = cellfun(Integration_fun, samples_cell, Projected_Distributions_cell, distribution_Sigma_cell, distribution_mu_c_cell, distribution_axes_cell, distribution_sigmae_cell, 'UniformOutput', false);
        
    %% Bias quantification %%
        if ED_Bias_Correction == true
            % The delta of the radius for which the expected distance is optimal for each distribution
            Bias_Plot = false;
            [prop_bias_cell, delta_radius_list] = Expected_Distance_Cylinder_Bias(distance_moment, cylinder_radius, Projected_Distributions_cell, alpha, number_samples, Bias_Plot, Diagnostics);
        else
            prop_bias_cell      = zeros(number_distributions, num_dim);
            delta_radius_list   = zeros(number_distributions, 1);
        end
        
        prop_bias_cell      = mat2cell(prop_bias_cell, ones(1, number_distributions), num_dim);
        delta_radius_cell   = num2cell(delta_radius_list);
        
    %% The expected distance %% 
        % Samples adjusted by the propagation bias
        Prop_Bias_fun   = @(samples_matrix, prop_bias) samples_matrix - prop_bias;
        samples_cell    = cellfun(Prop_Bias_fun, samples_cell, prop_bias_cell, 'UniformOutput', false);
    
        % Samples projected onto the cylinder axis
        Vector_Projection_fun       = @(samples_matrix) Point_to_Vector_Projection(samples_matrix, cylinder_direction, origin);     % Note that the distributions are centered w.r.t. the cylinder centre
        vector_proj_samples_cell    = cellfun(Vector_Projection_fun, samples_cell, 'UniformOutput', false);
    
        % The distances between each sample and the circle (centre at the origin)
        Distance_fun            = @(samples_matrix, proj_samples_matrix, delta_radius) abs(sqrt(sum((proj_samples_matrix - samples_matrix).^2, 2)) - (cylinder_radius + delta_radius)).^distance_moment;        % Note that the distance is taken to the given moment
        samples_distance_cell   = cellfun(Distance_fun, samples_cell, vector_proj_samples_cell, delta_radius_cell, 'UniformOutput', false);
        
        % The resulting expected distance for each ellipsoid
        Expected_Distance_fun   = @(Statistical_Outputs, samples_distance_list) Expected_Distance(Statistical_Outputs, samples_distance_list);
        expected_distance_list  = cellfun(Expected_Distance_fun, Statistical_Outputs_cell, samples_distance_cell);
       
        % And the mean value for all points
        expected_distance = mean(expected_distance_list);
        
    %% The dimensionally normalised expected distance %%
        % The expected distance is normalised by said value for a straight line, as curvature is unknown
        % As the rotation is similarly unknown, the expected value for any rotation is used which is approximately the average of the uncertainty of the 2D projection
        Line_Expected_Distance_fun      = @(Projected_Distributions) (sum(Projected_Distributions.Plane.Original.sigmae{1}) / sqrt(2*pi))^distance_moment;      % Note that it is taken to the desired moment
        line_expected_distance_list     = cellfun(Line_Expected_Distance_fun, Projected_Distributions_cell);
        
        norm_expected_distance_list     = expected_distance_list ./ line_expected_distance_list;
        
        % The mean value for the whole point cloud
        norm_expected_distance          = mean(norm_expected_distance_list);
        
    %% Result evaluation %%
        %--% Printed statements %--%
        if Print == true
            fprintf('E[D^%g] = %.3g m^%g \n', distance_moment, expected_distance, distance_moment);
            fprintf('E[D_N^%g] = %.3g \n', distance_moment, norm_expected_distance);
        end

        %--% Plot %--%
        if Plot == true
            % Number of coordinates used to discretise the cylinder and distributions
            number_coord = 1e2;

            % Estimate of the cylinder length
            distribution_mu_c_matrix    = vertcat(distribution_mu_c_cell{:});
            norm_list                   = sqrt(sum(distribution_mu_c_matrix.^2, 2));
            cylinder_length             = 2*max(norm_list);

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on

            % The cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.10, 'LineWidth', 2, 'DisplayName', 'Cylinder');

            % Its caps
            pl_top = plot3(cylinder_coord_x(1, :), cylinder_coord_y(1, :), cylinder_coord_z(1, :), 'LineWidth', 2, 'color', 'k');
            pl_bot = plot3(cylinder_coord_x(2, :), cylinder_coord_y(2, :), cylinder_coord_z(2, :), 'LineWidth', 2, 'color', 'k');            
            pl_top.HandleVisibility = 'Off';
            pl_bot.HandleVisibility = 'Off';
            
            % The distributions
            distribution_mu_matrix = vertcat(distribution_mu_cell{:});
            scatter3(distribution_mu_matrix(:, 1), distribution_mu_matrix(:, 2), distribution_mu_matrix(:, 3), 'filled', 'MarkerFaceColor', 'b', 'DisplayName', '\mu');

            for d = 1 : number_distributions
                % Distribution's properties
                distribution_mu     = distribution_mu_cell{d};
                distribution_sigmae = distribution_sigmae_cell{d};
                distribution_axes   = distribution_axes_cell{d};

                vector_proj_mu_t    = Projected_Distributions_cell{d}.Vector.Original.mu{1};
                vector_proj_mu      = vector_proj_mu_t + cylinder_centre;

                expected_distance   = expected_distance_list(d);
                expected_distance_1 = expected_distance^(1/distance_moment);        % s.t. its units is in metres

                % Surface
                [ellipsoid_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_sigmae, distribution_axes, number_coord);

                x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.25, 'DisplayName', sprintf('Point cloud, %s = 1', '\sigma'));

                % Vector from mu to its projection onto the cylinder axis, with direction to the cylinder surface
                orth_axis_vector    = vector_proj_mu - distribution_mu;
                orth_axis_vector_n  = orth_axis_vector / norm(orth_axis_vector);

                if norm(orth_axis_vector) < cylinder_radius
                    orth_axis_vector_n = -orth_axis_vector_n;
                end

                % Expected distance
                pl_ed = plot3(distribution_mu(1) + expected_distance_1 * [0, orth_axis_vector_n(1)], distribution_mu(2) + expected_distance_1 * [0, orth_axis_vector_n(2)], distribution_mu(3) + expected_distance_1 * [0, orth_axis_vector_n(3)], 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('E[D^%g]^{1/%g}', distance_moment, distance_moment));

                if d > 1
                    surf_el.HandleVisibility    = 'Off';
                    pl_ed.HandleVisibility      = 'Off';
                end
            end

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            axis equal
            view(45, 45);

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            % Legend
            legend('show', 'location', 'eastoutside');

            % Pause
            disp('The expected distance script will finish and the figure will close once a key is pressed');
            pause()

            close(1);
        end
    
    %% Expected distance function %%
        function E_D = Expected_Distance(Statistical_Outputs, samples_distance_list)
            % The probability-density weighted equation breaks if all samples have zero probability
            sample_density_list = Statistical_Outputs.Plane.probability_density;
            
            if sum(sample_density_list) == 0
                E_D = mean(samples_distance_list);                                                          % In that case, the average distance is used
            else
                E_D = sum(sample_density_list .* samples_distance_list) / sum(sample_density_list);        % Note that the finite confidence interval is taken into account
            end
        end    
end


