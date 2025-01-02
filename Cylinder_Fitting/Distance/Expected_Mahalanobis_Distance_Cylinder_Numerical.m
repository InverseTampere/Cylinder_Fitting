% The expected (squared) Mahalanobis distance over the given cylinder's cross-section is computed numerically using linear sampling over the circle
% Note that the bias can be taken into account by specifying the optional delta in radius

function [expected_Mahal_distance, expected_Mahal_dist_list] = Expected_Mahalanobis_Distance_Cylinder_Numerical(cylinder_centre, cylinder_radius, cylinder_direction, distance_moment, Point_Cloud_Distributions, Statistical_Values, Plot)
    
    %% Inputs %%
        % Point cloud distributions
        number_distributions    = Point_Cloud_Distributions.number_distributions;

        % Statistical values
        Confidence_interval     = Statistical_Values.Confidence_interval;

    %% Manual inputs %%
        % Projection onto the cross-section
        Proj_Diagnostics        = false;        % [true, false]

        % Determining the E[M^2] distance to a circle
        number_samples          = 1e1;          % [-] Integer of the number of samples used to numerically approximate it
        Circle_Sampling         = true;        % [true, false] To speed up computation, the distance can be estimated numerically by sampling the circle
        Circle_Diagnostics      = false;        % [true, false] Can be used to determine whether the estimate is correct or not

        % Correction of the radius bias due to curvature
        Radius_Bias_Correction  = false;        % [true, false] Whether or not radius bias is determined and corrected for
        
        Radius_Bias_Diagnostics = false;        % [true, false] Shows a plot of the optimal radius for each distribution
        GSS_Diagnostics         = false;        % [true, false] Can be used to determine whether GSS minimisation is working
        convergence_threshold   = 1e-2;         % [-] GSS convergence threshold. In terms of the squared Mahalanobis distance
        radius_factor           = 1.25;         % [-] Factor by which the radius is multipled/divided to denote the search bounds for GSS.
        max_iterations          = 1e1;          % [-] Maximum number of iterations for GSS.
               
    %% Projection onto the cross-section %%
        % The cross-sectional axes are unimportant
        [~, cylinder_vector_basis] = Rotation_3D(cylinder_centre, cylinder_direction, cylinder_centre);

        % Distributions projected on the cross-sectional plane   
        Projected_Distributions = Multivariate_Normal_Plane_Projection(cylinder_vector_basis, cylinder_centre, Point_Cloud_Distributions, Proj_Diagnostics);
        proj_distr_mu_cell      = Projected_Distributions.Plane.Projection.mu;
        proj_distr_sigmae_cell  = Projected_Distributions.Plane.Projection.sigmae;
        proj_distr_axes_cell    = Projected_Distributions.Plane.Projection.distr_axes;

        % The cylinder centre
        circle_centre_3D    = (cylinder_vector_basis' * cylinder_centre')';
        circle_centre       = circle_centre_3D(1 : 2);

    %% Quantifying the radius bias %%
        % The curvature of a circle results in an optimal radius being slightly outwards of the true one
        if Radius_Bias_Correction == true
            % The optimal radius is that for which the expected distance is minimal when the distribution lies on the circle
            Optimal_Radius_fun  = @(distr_mu, distr_axes, distr_sigmae) Optimal_Radius(cylinder_radius, circle_centre, distr_mu, distr_axes, distr_sigmae, radius_factor, convergence_threshold, max_iterations, GSS_Diagnostics);
            opt_radius_list     = cellfun(Optimal_Radius_fun, proj_distr_mu_cell, proj_distr_axes_cell, proj_distr_sigmae_cell);

            if Radius_Bias_Diagnostics == true
                for d = 1 : number_distributions
                    % This distribution's data
                    distr_mu        = proj_distr_mu_cell{d};
                    distr_axes      = proj_distr_axes_cell{d};
                    distr_sigmae    = proj_distr_sigmae_cell{d};
                    opt_radius      = opt_radius_list(d);

                    %--% Plot %--%
                    figure(d);
                    % Size and white background
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                    set(gcf, 'color', [1, 1, 1])     
        
                    hold on
                    grid on

                    % True radius
                    number_coord                    = 1e2;
                    [x_circle_list, y_circle_list]  = Circle_Coordinates(cylinder_radius, 0, 0, number_coord);
                    plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'k', 'DisplayName', 'True radius');

                    % Optimal radius circle
                    [x_circle_list, y_circle_list]  = Circle_Coordinates(opt_radius, 0, 0, number_coord);
                    plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('Optimal radius = %.3g m', opt_radius));

                    % Distribution lying on the circle
                    circle_distr_mu     = (distr_mu - circle_centre) / norm(distr_mu - circle_centre) * cylinder_radius;
                    distr_coord_matrix  = Ellipse_Coordinate_Generator(circle_distr_mu, distr_axes, distr_sigmae, number_coord);

                    plot(distr_coord_matrix(:, 1), distr_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', sprintf('Distribution (1%s)', '\sigma'));

                    % Axes
                    xlabel('x [m]');
                    ylabel('y [m]');
                    axis equal

                    % Legend
                    legend('show', 'location', 'eastoutside');

                    set(gca, 'FontSize', 15);
                    set(gca, 'LineWidth', 2);
                end

                % Pause message
                disp('The optimal radius has been determined. The figures will close and script continue upon a key-press.');
                pause();

                close(1:d);
            end
        else
            % If bias correction does not take place, the optimal radius is taken as the true one
            opt_radius_list = repmat(cylinder_radius, [number_distributions, 1]);
        end

        % Converted to cell array for cellfun
        opt_radius_cell = num2cell(opt_radius_list);

    %% The expected Mahalanobis distance %%
        % The expected Mahalanobis distance
        Circle_Mahal_Distance_fun   = @(distr_mu, distr_axes, distr_sigmae, circle_radius) Expected_Mahalanobis_Distance_Circle_Numerical(circle_centre, circle_radius, distr_mu, distr_sigmae, distr_axes, distance_moment, number_samples, Confidence_interval, Circle_Sampling, Circle_Diagnostics);
        expected_Mahal_dist_list    = cellfun(Circle_Mahal_Distance_fun, proj_distr_mu_cell, proj_distr_axes_cell, proj_distr_sigmae_cell, opt_radius_cell);

        % The expected distance as a whole is taken as the average
        expected_Mahal_distance = mean(expected_Mahal_dist_list);

    %% Plot %%
        if Plot == true
            figure(1);
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on

            % Circle
            number_coord                    = 1e2;
            [x_circle_list, y_circle_list]  = Circle_Coordinates(cylinder_radius, circle_centre(1), circle_centre(2), number_coord);
            plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Circle radius');

            % Distributions
            for d = 1 : number_distributions
                % This distribution's info
                distr_mu        = proj_distr_mu_cell{d};
                distr_axes      = proj_distr_axes_cell{d};
                distr_sigmae    = proj_distr_sigmae_cell{d};
                Mahal_dist      = expected_Mahal_dist_list(d);

                % Distribution at 1 sigma
                distr_coord_matrix  = Ellipse_Coordinate_Generator(distr_mu, distr_axes, distr_sigmae, number_coord);
                pl_distr            = plot(distr_coord_matrix(:, 1), distr_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', sprintf('Distribution (1%s)', '\sigma'));

                % scaled Mahalanobis distance
                scaled_Mahal_dist   = sqrt(prod(distr_sigmae)) * Mahal_dist^(1/distance_moment);        
                vector_sign         = sign(cylinder_radius - norm(distr_mu - circle_centre));           % Negative if the vector should point inwards
                distance_vector     = scaled_Mahal_dist * vector_sign * (distr_mu - circle_centre) / norm(distr_mu - circle_centre);

                pl_dist = plot(distr_mu(1) + [0, distance_vector(1)], distr_mu(2) + [0, distance_vector(2)], 'LineWidth', 2, 'color', 'k', 'DisplayName', sprintf('(|%s|E[M^2])^{1/2}', '\Sigma'));

                if d > 1
                    pl_distr.HandleVisibility   = 'Off';
                    pl_dist.HandleVisibility    = 'Off';
                end
            end

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            axis equal

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            % Legend
            legend('show', 'location', 'eastoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            % Pause message
            fprintf('E[M^%i] = %.3g. The figure will close and script end upon a key-press. \n', distance_moment, expected_Mahal_distance);
            pause();

            close(1);
        end

    %% Local functions %%
        % The optimal radius of a given distribution
        function opt_radius = Optimal_Radius(cylinder_radius, circle_centre, distr_mu, distr_axes, distr_sigmae, radius_factor, convergence_threshold, max_iterations, GSS_Diagnostics)
            % Golden-section search for minimisation
            radius_lower_bound = cylinder_radius / radius_factor;
            radius_upper_bound = cylinder_radius * radius_factor;

            Circle_Distance_fun = @(radius) Expected_Mahal_Distance_Circle_Radius(radius, cylinder_radius, circle_centre, distr_mu, distr_axes, distr_sigmae);
            opt_radius          = Golden_Section_Search(Circle_Distance_fun, radius_lower_bound, radius_upper_bound, convergence_threshold, max_iterations, GSS_Diagnostics, GSS_Diagnostics);
        end
        
        % The squared expected Mahalanobis distance of a distribution lying on the true circle for a circle of a different radius
        function expected_Mahal_distance = Expected_Mahal_Distance_Circle_Radius(fake_radius, true_radius, circle_centre, distr_mu, distr_axes, distr_sigmae)
            % The distribution's centre is displaced to lie on the true circle, which is placed at the origin
            circle_mu   = (distr_mu - circle_centre) / norm(distr_mu - circle_centre) * true_radius;

            num_dim     = length(circle_mu);
            origin      = zeros(1, num_dim);

            % E[M^m] of the distribution lying on the true circle
            expected_Mahal_distance = Expected_Mahalanobis_Distance_Circle_Numerical(origin, fake_radius, circle_mu, distr_sigmae, distr_axes, distance_moment, number_samples, Confidence_interval, Circle_Sampling, Circle_Diagnostics);
        end
end