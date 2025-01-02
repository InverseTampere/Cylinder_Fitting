% An infinite plane is fitted to the given point cloud distribution by minimising the expected squared Mahalanobis distance. 
% The normal vector may be fixed or optimised for, in which case the given vector is used as first estimate
% The height is with respect to the origin along the normal vector, resulting in the plane point
% The given bounds margin is a factor of pi/2 for the elevation angle if the normal vector is optimised for 

function [Plane_Geometry, plane_normal_vector, plane_point, plane_height, expected_Mahal_dist_plane, Optimiser_Diagnostics] = Fuzzy_Infinite_Plane_Fitting(Point_Cloud_Distributions, Fitting_Parameters, Fixed_normal_vector, plane_normal_vector, Print, Plot)

    %% Structure inputs %%
        % Fitting parameters
        bounds_margin           = Fitting_Parameters.bounds_margin;
        point_weights_list      = Fitting_Parameters.point_weights_list;

    %% Manual inputs %%
        % Optimiser options for determination of the normal vector
        max_iterations          = 1e2;                  % [-] Maximum number of iterations by the optimiser. If it is reached, a warning is displayed
        function_tolerance      = 1e-4;                 % [-] Minimum size for the step in the objective function. Note that it is normalised by the initial value
        step_tolerance          = 1e-4;                 % [-] Minimum size for the step in the design parameters. Note that they are normalised between [0, 1]
        FD_step_size            = 1e-4;                 % [-] Step size taken to evaluate the finite difference. Note that sqrt(eps) is the default and that it doesn't affect the analytical gradient or Hessian
        
        % Plot settings
        Proj_Diagnostics        = false;                % [true, false] For diagnostic purposes projection onto the normal vector can be shown
        Mahal_Diagnostics       = false;                % [true, false] Shows a plot to diagnose the expected squared Mahalanobis distance with

    %% Optimal plane height for given normal vector %%
        % Spherical angles of the given normal vector
        [~, azim_angle, elev_angle, ~, ~, ~]    = Vector_Spherical_Angle_Conversion(plane_normal_vector, [], []);
        [azim_ind, elev_ind]                    = deal(1, 2);
        Parameter_Indices                       = struct('azim_angle', azim_ind, 'elev_angle', elev_ind);

        spherical_angles                        = zeros(1, 2);
        spherical_angles([azim_ind, elev_ind])  = [azim_angle, elev_angle];

        % Initial expected Mahalanobis distance
        [expected_Mahal_dist_plane, plane_height] = Plane_Expected_Mahal_Distance(spherical_angles, point_weights_list, Point_Cloud_Distributions, Proj_Diagnostics, Mahal_Diagnostics);

    %% Optimised normal vector %%
        if Fixed_normal_vector == false
            %--% Optimiser settings %--%
            % Optimisation bounds
            elev_deviation      = bounds_margin * pi/2;
            elev_LB             = max(-pi/2, elev_angle - elev_deviation);
            elev_UB             = min(pi, elev_angle + elev_deviation);
            [azim_LB, azim_UB]  = deal(-pi, pi);
           
            % The combined sets of bounds
            spherical_LB                        = zeros(1, 2);
            spherical_LB([azim_ind, elev_ind])  = [azim_LB, elev_LB];
            spherical_UB                        = zeros(1, 2);
            spherical_UB([azim_ind, elev_ind])  = [azim_UB, elev_UB];
        
            % The optimisation process is saved in matrices
            [spherical_angle_steps, Gradient_steps, Objective_value_steps] = deal([]);

            % Optimiser options
            Options = optimoptions('fmincon', 'OutputFcn', @Optim_Path, 'Algorithm', 'interior-point', 'SubProblemAlgorithm', 'cg', 'EnableFeasibilityMode', true, 'MaxIterations', max_iterations, 'FunctionTolerance', function_tolerance, 'StepTolerance', step_tolerance, 'FiniteDifferenceStepSize', FD_step_size, 'Display', 'off');

            %--% Optimal plane normal vector %--%
            % Optimal spherical angles
            Plane_Fitting_Objective_fun = @(spherical_angles) Plane_Expected_Mahal_Distance(spherical_angles, point_weights_list, Point_Cloud_Distributions, Proj_Diagnostics, Mahal_Diagnostics) / expected_Mahal_dist_plane;
            
            [opt_spherical_angles, ~, ~, ~, ~, Optimum_gradient_list, Optimum_Hessian_matrix] = fmincon(Plane_Fitting_Objective_fun, spherical_angles, [], [], [], [], spherical_LB, spherical_UB, [], Options);
            
            % Conversion to vector
            [azim_angle, elev_angle]                = deal(opt_spherical_angles(azim_ind), opt_spherical_angles(elev_ind));
            [plane_normal_vector, ~, ~, ~, ~, ~]    = Vector_Spherical_Angle_Conversion([], azim_angle, elev_angle);

            % Corresponding plane height
            [expected_Mahal_dist_plane, plane_height] = Plane_Expected_Mahal_Distance(opt_spherical_angles, point_weights_list, Point_Cloud_Distributions, Proj_Diagnostics, Mahal_Diagnostics);

            %--% Checks %--%
            % Check if the maximum number of iterations was reached
            num_optimiser_steps = length(Objective_value_steps);
    
            if num_optimiser_steps == max_iterations
                warning('The maximum number of iterations for the optimiser was reached. Check whether this should be increased or whether the objective function is poor.')
            end
    
            % Check if the objective function is convex at the solution (i.e. if all eigenvalues of the Hessian are positive)
            % Note that as the optimising algorithm is quasi-Newton, convexity is not guaranteed
            Hessian_eigenvalues = eig(Optimum_Hessian_matrix);
            
            if min(Hessian_eigenvalues) < 0
                fprintf('The objective function is not convex at the found solution. The lowest Hessian eigenvalue is %.3g \n', min(Hessian_eigenvalues));
            end

            %--% Optimiser steps %--%
            % Step sizes 
            Parameter_steps_delta_matrix    = diff(spherical_angle_steps, 1, 1);
            Optimiser_step_size_list        = [0; sqrt(sum(Parameter_steps_delta_matrix.^2, 2))];             % A 0 value is appended for the first step
            
            % Labels and units for each parameter
            parameter_labels                        = cell(1, 2);
            parameter_labels([azim_ind, elev_ind])  = {'azim_angle', 'elev_angle'};
    
            parameter_units                         = {'rad', 'rad'};
    
            % Structures for the (gradient) steps
            Geometry_Steps = struct();
            Gradient_Steps = struct();
    
            for i = 1 : 2
                parameter_label                     = parameter_labels{i};
                Geometry_Steps.(parameter_label)    = spherical_angle_steps(:, i);
                Gradient_Steps.(parameter_label)    = Gradient_steps(:, i);
            end
    
            % Complete structures
            Optimum_First_Derivatives   = struct('first_derivative_list', Optimum_gradient_list');              % Note that they are changed to a row vector
            Optimum_Second_Derivatives  = struct('second_derivative_list', diag(Optimum_Hessian_matrix)');
            Optimum_Hessian             = struct('Hessian_matrix', Optimum_Hessian_matrix);
    
            Optimiser_Steps             = struct('Geometry_Steps', Geometry_Steps, 'Gradient_Steps', Gradient_Steps, 'Step_sizes', Optimiser_step_size_list, 'Objective_value_steps', Objective_value_steps, 'num_optimiser_steps', num_optimiser_steps);
            Optimiser_Diagnostics       = struct('parameter_labels', {parameter_labels}, 'parameter_units', {parameter_units}, 'num_parameters', length(parameter_units), 'Parameter_Indices', Parameter_Indices, 'Optimum_First_Derivatives', Optimum_First_Derivatives, 'Optimum_Second_Derivatives', Optimum_Second_Derivatives, 'Optimum_Hessian', Optimum_Hessian, 'Optimiser_Steps', Optimiser_Steps);

        else
            Optimiser_Diagnostics       = [];           % Otherwise no optimisation is performed
        end

        % A point lying on the plane
        plane_point = plane_normal_vector * plane_height;

        % For an infinite plane the orientation of the planar vectors is not so important
        num_dim                     = length(plane_normal_vector);
        origin                      = zeros(1, num_dim);
        [~, plane_vector_basis, ~]  = Vector_Based_Rotation(origin, plane_normal_vector, origin);
        
        % Geometry structure
        Plane_Geometry = struct('plane_point', plane_point, 'plane_height', plane_height, 'plane_normal_vector', plane_normal_vector, 'plane_vector_basis', plane_vector_basis, 'azim_angle', azim_angle, 'elev_angle', elev_angle);

    %% Fitting information %%
        % Printed messages
        if Print == true
            disp('The plane has been fitted with:')
            fprintf('h = %.3g, n_vec = [%.3g, %.3g, %.3g], E[M^2] = %.3g \n', plane_height, plane_normal_vector, expected_Mahal_dist_plane)
        end

        % Plot
        if Plot == true
            % Plot settings
            m_STD           = 1;
            number_coord    = 1e2;

            % Distribution properties
            number_distributions        = Point_Cloud_Distributions.number_distributions;
            distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;
            distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
            distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
            point_cloud_matrix          = vertcat(distribution_mu_cell{:});

            %--% Plot %--%
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])

            hold on
            grid on            

            % Expected values of the distributions
            scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 50, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', sprintf('%s', '\mu'));

            for i = 1 : number_distributions
                % Distribution properties
                distribution_mu     = distribution_mu_cell{i};
                distribution_sigmae = distribution_sigmae_cell{i};
                distribution_radii  = distribution_sigmae * m_STD;
                distribution_axes   = distribution_axes_cell{i};

                % Its coordinates at m sigma
                [distr_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_coord);

                x_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                y_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                z_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                % Its surface patch
                surf_distr_name = sprintf('Distribution, %.2g %s', m_STD, '\sigma');
                distr_surf = surf(x_distribution, y_distribution, z_distribution, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.25, 'DisplayName', surf_distr_name);

                if i > 1
                    distr_surf.HandleVisibility = 'Off';
                end
            end

            % The fitted plane
            plane_corner_matrix = Plane_Corner_Points(plane_normal_vector, plane_point, point_cloud_matrix);
            patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'm', 'FaceAlpha', 0.25, 'DisplayName', 'Fitted plane');   

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            axis equal

            % Viewing angle
            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Pause message
            disp('The plane has been fitted. The figure will close and script will end upon a key-press.');
            pause();

            close(1);
        end

    %% Local functions %%
        % E[M^2] of the given plane
        function [expected_Mahal_dist_plane, plane_height] = Plane_Expected_Mahal_Distance(spherical_angles, point_weights_list, Point_Cloud_Distributions, Proj_Diagnostics, Mahal_Diagnostics)
            % Spherical angles
            azim_angle = spherical_angles(azim_ind);
            elev_angle = spherical_angles(elev_ind);
        
            % Conversion to plane normal vector 
            [plane_normal_vector, ~, ~, ~, ~, ~] = Vector_Spherical_Angle_Conversion([], azim_angle, elev_angle);

            % Projection vector basis
            num_dim                     = length(plane_normal_vector);
            origin                      = zeros(1, num_dim);
            [~, plane_vector_basis, ~]  = Vector_Based_Rotation(origin, plane_normal_vector, origin);

            % Projection onto the normal vector
            Plane_Proj_Distributions    = Multivariate_Normal_Plane_Projection(plane_vector_basis, origin, Point_Cloud_Distributions, Proj_Diagnostics);
            proj_mu_list                = vertcat(Plane_Proj_Distributions.Vector.Projection.mu{:});
            proj_sigma_list             = vertcat(Plane_Proj_Distributions.Vector.Projection.sigma{:});

            % Optimal plane height
            plane_height = sum(proj_mu_list ./ proj_sigma_list.^2) / sum(proj_sigma_list.^(-2));

            % Expected Mahalanobis distance
            expected_Mahal_dist_plane_list      = proj_sigma_list.^(-2)*plane_height.^2 - 2*proj_mu_list.*proj_sigma_list.^(-2)*plane_height + proj_mu_list.^2.*proj_sigma_list.^(-2) + 1;
            expected_Mahal_dist_plane_list_w    = point_weights_list .* expected_Mahal_dist_plane_list; 
            expected_Mahal_dist_plane           = mean(expected_Mahal_dist_plane_list_w);
    
            % Diagnostics plot 
            if Mahal_Diagnostics == true   
                % Plot settings
                m_STD           = 1;
                number_coord    = 1e2;

                % Distribution properties
                number_distributions        = Point_Cloud_Distributions.number_distributions;
                distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;
                distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
                distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
                point_cloud_matrix          = vertcat(distribution_mu_cell{:});
    
                % Rescaled Mahalanobis distance
                scaled_expected_Mahal_dist_list = sqrt(expected_Mahal_dist_plane_list) .* proj_sigma_list;
    
                %--% Plot %--%
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])
    
                hold on
                grid on
    
                % The expected values
                scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'DisplayName', sprintf('%s', '\mu'));
    
                for d = 1 : number_distributions
                    % Distribution properties
                    distribution_mu     = distribution_mu_cell{d};
                    distribution_sigmae = distribution_sigmae_cell{d};
                    distribution_radii  = distribution_sigmae * m_STD;
                    distribution_axes   = distribution_axes_cell{d};
                    scaled_Mahal_dist   = scaled_expected_Mahal_dist_list(d);
    
                    % Its coordinates at m sigma
                    [distr_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_coord);
    
                    x_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                    y_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                    z_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);
    
                    % Its surface patch
                    surf_distr_name = sprintf('Distribution, %.2g %s', m_STD, '\sigma');
                    distr_surf = surf(x_distribution, y_distribution, z_distribution, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', surf_distr_name);
    
                    % The scaled expected Mahalanobis distance
                    dist_vector = sign(plane_height - proj_mu_list(d)) * plane_normal_vector * scaled_Mahal_dist;
                    dist_plot   = plot3(distribution_mu(1) + [0, dist_vector(1)], distribution_mu(2) + [0, dist_vector(2)], distribution_mu(3) + [0, dist_vector(3)], 'LineWidth', 2, 'color', 'k', 'DisplayName', sprintf('sqrt(E[M^2]) * %s', '\sigma_h'));
    
                    if d > 1
                        distr_surf.HandleVisibility = 'Off';
                        dist_plot.HandleVisibility  = 'Off';
                    end
                end
    
                % The fitted plane
                plane_point         = plane_height * plane_normal_vector;                                               % Note that the plane height is w.r.t. the origin
                plane_corner_matrix = Plane_Corner_Points(plane_normal_vector, plane_point, point_cloud_matrix);
                patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'm', 'FaceAlpha', 0.25, 'DisplayName', 'Fitted plane');   
    
                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');
    
                axis equal
    
                % Viewing angle
                view(45, 45);
    
                % Legend
                legend('show', 'location', 'northoutside');
    
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off   
    
                % Pause message
                disp('The expected Mahalanobis distance to the plane has been computed. The plot will close and script end upon a key-press.');
                pause();
    
                close(1);
            end
        end

        % Optimisation path function
        function stop = Optim_Path(spherical_angles, Optimisation_Values, ~)
            stop = false;

            % Concatenate objective value, geometry and various other terms computed during each iteration
            spherical_angle_steps   = [spherical_angle_steps; spherical_angles];
            Gradient_steps          = [Gradient_steps; Optimisation_Values.gradient'];
            Objective_value_steps   = [Objective_value_steps; Optimisation_Values.fval];
        end
end
