% This function determines the cylinder cross-section and direction through minimising the expected Mahalanobis distance using fmincon

function [Optimal_Geometry, Init_Geometry, Objective_value, Point_Cloud_Distributions, Optimiser_Diagnostics] = Cylinder_Cross_Section_and_Direction_Estimation_fmincon(Point_Cloud_Coord, Scanning_Parameters, Scanner_Parameters, Statistical_Values, Fitting_Parameters, Output_Decisions)

    %% Inputs %% 
        % Scanning parameters
        Scanner_loc_cell            = Scanning_Parameters.Scanner_loc_cell;
        number_scanners             = Scanning_Parameters.number_scanners;

        % Scanner parameters
        range_bias                  = Scanner_Parameters.range_bias;     
        
        % Fitting parameters
        RANSAC                      = Fitting_Parameters.RANSAC;
        bounds_margin               = Fitting_Parameters.bounds_margin;
        max_UU_iterations           = Fitting_Parameters.max_UU_iterations;
        UU_convergence_threshold    = Fitting_Parameters.UU_convergence_threshold;

        % Statistical values
        Confidence_interval         = Statistical_Values.Confidence_interval;
        
        % Output decisions
        Print                       = Output_Decisions.Print;
        Plot                        = Output_Decisions.Plot;

    %% Manual inputs %%
        % If desired, intermediate results can be checked
        LS_Plot                 = false;        % [true, false]
        LS_Diagnostics          = false;        % [true, false]
        Transform_Diagnostics   = false;        % [true, false]
        Distr_Diagnostics       = false;        % [true, false]
        Suppress_Warnings       = true;         % [true, false] Custom warnings may be suppressed if desired

        % Optimiser inputs
        max_iterations          = 1e2;          % [-] Maximum number of iterations by the optimiser. If it is reached, a warning is displayed
        function_tolerance      = 1e-4;         % [-] Minimum size for the step in the objective function. Note that it is normalised by the initial value
        step_tolerance          = 1e-4;         % [-] Minimum size for the step in the design parameters. Note that they are normalised between [0, 1]
        FD_step_size            = 1e-4;         % [-] Step size taken to evaluate the finite difference. Note that sqrt(eps) is the default and that it doesn't affect the analytical gradient or Hessian

    %% Initial geometry estimate %%
        % Least-squares is used for the initial geometry estimate, either through RANSAC or not
        if RANSAC == true
            Init_Geometry = RANSAC_Cylinder_Fitting(Point_Cloud_Coord, Scanning_Parameters, Scanner_Parameters, Fitting_Parameters, Statistical_Values, LS_Plot);
        else
            Init_Geometry = Least_Squares_Cylinder_Fitting(Point_Cloud_Coord, Scanning_Parameters, Fitting_Parameters, Confidence_interval, LS_Plot, LS_Diagnostics);
        end
        
        [cylinder_centre_init, cylinder_direction_init, cylinder_radius_init, cylinder_length] = deal(Init_Geometry.Cylinder_centre, Init_Geometry.Cylinder_direction, Init_Geometry.Cylinder_radius, Init_Geometry.Cylinder_length);

        % The third component is made to be positive by flipping the vector if need be, as sin(elev) >= 0
        num_dim                 = length(cylinder_direction_init);
        cylinder_direction_init = sign(cylinder_direction_init(num_dim)) * cylinder_direction_init;

    %% Translation and rotation for robustness %%
        % To increase robustness, everything is centered on the centroid of the given point cloud
        [Point_Cloud_Coord_c, ~, ~, point_cloud_centroid, Scanning_Parameters_c, ~, ~] = Point_Cloud_Centering(Point_Cloud_Coord, [], Scanning_Parameters);

        % Additionally, the data is rotated s.t. the elevation angle of the cylinder axis is pi/4 to avoid gimbal lock and such that the initial elevation angle is in the middle of the optimisation range
        origin                          = zeros(1, num_dim);
        [~, cylinder_vector_basis, ~]   = Vector_Based_Rotation(origin, cylinder_direction_init, origin);
        
        pitch                           = pi/4;
        pitched_vector_basis            = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];

        gimbal_rotation_matrix          = pitched_vector_basis * cylinder_vector_basis;

        [Point_Cloud_Coord_r, ~, ~, Scanning_Parameters_r, Scanner_loc_cell_r, ~] = Point_Cloud_Rotation(gimbal_rotation_matrix, Point_Cloud_Coord_c, [], Scanning_Parameters_c);

        % The initial geometry must be translated and rotated as well
        cylinder_direction_init_r               = (gimbal_rotation_matrix * cylinder_direction_init')';
        cylinder_centre_init_c                  = cylinder_centre_init - point_cloud_centroid;
        cylinder_centre_init_r                  = (gimbal_rotation_matrix * cylinder_centre_init_c')';
        [circle_centre_init_matrix_r, ~, ~, ~]  = Circle_Cylinder_Centre_Conversion([], [], cylinder_centre_init_r, cylinder_direction_init_r, origin, Scanner_loc_cell_r);
        circle_centre_init_r                    = circle_centre_init_matrix_r(1, :);

        % Diagnostic plot to check the transformation
        if Transform_Diagnostics == true
            % The number of coordinates in each ellipsoid and the cylinder
            number_coord = 1e3;     
            
            scanner_colours = cbrewer('qual', 'Set2', max(number_scanners, 3));     % The colours used for the scanners
            scanner_colours = max(scanner_colours, 0);
            scanner_colours = min(scanner_colours, 1);
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            %--% Original coordinate frame %--%
            subplot(1, 2, 1)
            hold on
            grid on
                
            % The initial cylinder surface
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius_init, cylinder_length, cylinder_centre_init, cylinder_direction_init, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial cylinder');
  
            % The point cloud
            point_cloud_cell = Point_Cloud_Coord.point_cloud_cell;

            for s = 1 : number_scanners
                % This scanner's point cloud
                scanner_colour      = scanner_colours(s, :);
                point_cloud_matrix  = point_cloud_cell{s};

                scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 2), 'MarkerFaceColor', scanner_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Point cloud');
            end

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')
            axis equal
            
            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            %--% Transformed coordinate frame %--%
            subplot(1, 2, 2)
            hold on
            grid on
                
            % The initial cylinder surface
            [init_cylinder_coord_a, init_cylinder_coord_b, init_cylinder_coord_c, ~] = Cylinder_Surface_Generator(cylinder_radius_init, cylinder_length, cylinder_centre_init_r, cylinder_direction_init_r, number_coord);
            surf(init_cylinder_coord_a, init_cylinder_coord_b, init_cylinder_coord_c, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial cylinder');
  
            % The point cloud
            point_cloud_cell = Point_Cloud_Coord_r.point_cloud_cell;

            for s = 1 : number_scanners
                % This scanner's point cloud
                scanner_colour      = scanner_colours(s, :);
                point_cloud_matrix  = point_cloud_cell{s};
                
                scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 2), 'MarkerFaceColor', scanner_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Point cloud');
            end

            % Axes
            xlabel('a [m]')
            ylabel('b [m]')
            zlabel('c [m]')
            axis equal
            
            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off          
        end

    %% Solution bounds %%
        % The geometry parameters' indices (as the optimiser's inputs must be vector format)
        centre_ind                  = 1 : 2;
        azim_ind                    = max(centre_ind) + 1;
        elev_ind                    = azim_ind + 1;
        radius_ind                  = elev_ind + 1;
        num_optim_geo_parameters    = radius_ind;

        Geometry_Indices            = struct('centre', centre_ind, 'azim_angle', azim_ind, 'elev_angle', elev_ind, 'radius', radius_ind);

        % Finally, the set of initial centered geometry parameters
        inf_cylinder_geometry_init                                                  = zeros(1, num_optim_geo_parameters);
        [azim_angle_init_r, elev_angle_init_r]                                      = deal(0, pitch);
        inf_cylinder_geometry_init([centre_ind, azim_ind, elev_ind, radius_ind])    = [circle_centre_init_r, azim_angle_init_r, elev_angle_init_r, cylinder_radius_init];

        % Reasonable bounds for the optimisation process are based on the least-squares geometry s.t. the normalised variables can be directly compared
        circle_centre_LB = circle_centre_init_r - bounds_margin * cylinder_radius_init;
        circle_centre_UB = circle_centre_init_r + bounds_margin * cylinder_radius_init;
        
        radius_LB = cylinder_radius_init / (1 + bounds_margin);
        radius_UB = cylinder_radius_init * (1 + bounds_margin);

        % The direction vectors can deviate by the radius times the margin, at the cylinder length's distance
        % Note that this deviation only results in allowed changes to the elevation angle. The azimuth angle remains unbounded
        elev_deviation      = atan(cylinder_radius_init * bounds_margin / cylinder_length);

        elev_LB             = max(-pi/2, elev_angle_init_r - elev_deviation);
        elev_UB             = min(pi, elev_angle_init_r + elev_deviation);
        [azim_LB, azim_UB]  = deal(-3*pi, 3*pi);

        % The combined sets of bounds
        geometry_LB                                                 = zeros(1, num_optim_geo_parameters);
        geometry_LB([centre_ind, azim_ind, elev_ind, radius_ind])   = [circle_centre_LB, azim_LB, elev_LB, radius_LB];

        geometry_UB                                                 = zeros(1, num_optim_geo_parameters);
        geometry_UB([centre_ind, azim_ind, elev_ind, radius_ind])   = [circle_centre_UB, azim_UB, elev_UB, radius_UB];

        % As a result, the optimisation bounds are between 0 and 1
        Optim_LB = zeros(1, num_optim_geo_parameters);
        Optim_UB = ones(1, num_optim_geo_parameters);

    %% Iterative shape fitting %%
        % Turning custom warnings off to avoid spam with updated uncertainty if desired
        if Suppress_Warnings == true && max_UU_iterations > 1
            warning('off', 'optimiser:max_iter');
            warning('off', 'optimiser:convexity');
            warning('off', 'bounds:centre');
            warning('off', 'bounds:radius');
            warning('off', 'bounds:axis');
        end

        % Certain information is kept through the iterations for diagnostic purposes
        [Infinite_cylinder_geometry_steps_cell_n, Gradient_steps_cell_n, Objective_value_steps_cell_n] = deal(cell(1, max_UU_iterations));

        inf_cylinder_geometry_n_matrix  = NaN(max_UU_iterations, num_optim_geo_parameters);
        avg_Sigma_det_list              = NaN(1, max_UU_iterations);
        number_optim_steps_list         = NaN(1, max_UU_iterations);
        geometry_difference_list        = NaN(1, max_UU_iterations);

        % Initialising the loop
        inf_cylinder_geometry_prev_n    = (inf_cylinder_geometry_init - geometry_LB) ./ (geometry_UB - geometry_LB);      % Normalised LS estimate
        Convergence                     = false;
        iter                            = 0;

        % Each iteration the distributions are updated based on the previous geometry
        while iter < max_UU_iterations && Convergence == false
            % Update the iteration counter
            iter = iter + 1;

            % Previous geometry estimate
            [~, cylinder_direction_r, ~, cylinder_centre_r, ~] = Geometry_Retrieval(inf_cylinder_geometry_prev_n, geometry_LB, geometry_UB, gimbal_rotation_matrix, point_cloud_centroid, Scanner_loc_cell_r, Geometry_Indices);

            %--% Uncertainty %--%
            % The range and radial uncertainty of each point
            [sigma_radial_cell, sigma_prop_cell, ~, ~] = Cylindrical_Object_Uncertainty(cylinder_centre_r, cylinder_direction_r, Scanner_loc_cell_r, Scanner_Parameters, Point_Cloud_Coord_r);

            % The 3D distributions
            alpha                       = 1 - Confidence_interval/100;          % The confidence interval is changed to alpha
            Point_Cloud_Distributions_r = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord_r, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell_r, range_bias, Distr_Diagnostics);

            distr_Sigma_cell            = Point_Cloud_Distributions_r.distribution_Sigma_cell;
            distr_Sigma_det_list        = cellfun(@det, distr_Sigma_cell);
            avg_Sigma_det_list(iter)    = mean(distr_Sigma_det_list);

            %--% Optimisation %--%
            % Centered geometry and a centered point cloud are used with the initial objective value for normalisation
            Objective_value_init = Cylinder_Cross_Section_Objective(inf_cylinder_geometry_prev_n, Geometry_Indices, geometry_LB, geometry_UB, 1, Fitting_Parameters, Point_Cloud_Distributions_r, Scanning_Parameters_r, Statistical_Values);

            Infinite_Cylinder_Optimisation_fun = @(inf_cylinder_geometry_n) Cylinder_Cross_Section_Objective(inf_cylinder_geometry_n, Geometry_Indices, geometry_LB, geometry_UB, Objective_value_init, Fitting_Parameters, Point_Cloud_Distributions_r, Scanning_Parameters_r, Statistical_Values);      

            % Optimisation with the interior-point algorithm
            Options = optimoptions('fmincon', 'OutputFcn', @Optim_Path, 'Algorithm', 'interior-point', 'SubProblemAlgorithm', 'cg', 'EnableFeasibilityMode', true, 'MaxIterations', max_iterations, 'FunctionTolerance', function_tolerance, 'StepTolerance', step_tolerance, 'FiniteDifferenceStepSize', FD_step_size, 'Display', 'off');

            [inf_cylinder_geometry_n, Objective_value, ~, ~, ~, Optimum_gradient_list_n, Optimum_Hessian_matrix_n] = fmincon(Infinite_Cylinder_Optimisation_fun, inf_cylinder_geometry_prev_n, [], [], [], [], Optim_LB, Optim_UB, [], Options);

            %--% Checks %--%
            % Check if the optimum is worse than the initial geometry, in which case sqp is used as the optimisation algorithm
            if Objective_value > 1
                % The previous optimisation run is removed from the history
                [Infinite_cylinder_geometry_steps_cell_n{iter}, Gradient_steps_cell_n{iter}, Objective_value_steps_cell_n{iter}] = deal([]);

                % Updated optimisation method
                Options = optimoptions('fmincon', 'OutputFcn', @Optim_Path, 'Algorithm', 'sqp', 'MaxIterations', max_iterations, 'FunctionTolerance', function_tolerance, 'StepTolerance', step_tolerance, 'FiniteDifferenceStepSize', FD_step_size, 'Display', 'off');                        % Algorithm is changed to sqp

                % New optimisation run
                [inf_cylinder_geometry_n, Objective_value, ~, ~, ~, Optimum_gradient_list_n, Optimum_Hessian_matrix_n] = fmincon(Infinite_Cylinder_Optimisation_fun, inf_cylinder_geometry_prev_n, [], [], [], [], Optim_LB, Optim_UB, [], Options);    % Optimisation with sqp
            end

            % The geometry is added to the matrix
            inf_cylinder_geometry_n_matrix(iter, :) = inf_cylinder_geometry_n;

            % Check if the maximum number of iterations was reached
            number_optimiser_steps          = length(Objective_value_steps_cell_n{iter});
            number_optim_steps_list(iter)   = number_optimiser_steps;

            if number_optimiser_steps == max_iterations
                warning('optimiser:max_iter', 'The maximum number of iterations for the optimiser was reached. Check whether this should be increased or whether the objective function is poor.')
            end

            % Check if any of the bounds are reached, which may indicate that the optimum lies outside the optimisation space
            [centre_n, azim_n, elev_n, radius_n] = deal(inf_cylinder_geometry_n(centre_ind), inf_cylinder_geometry_n(azim_ind), inf_cylinder_geometry_n(elev_ind), inf_cylinder_geometry_n(radius_ind));

            if max(centre_n) > 1 - 2*step_tolerance || min(centre_n) < 2*step_tolerance                                                                             % Note that the finite step size means the bounds are unlikely to be hit exactly
                warning('bounds:centre', 'The bounds for optimisation have restricted the optimal centre location. Try increasing the bounds margin.');
            end
            if radius_n > 1 - 2*step_tolerance || radius_n < 2*step_tolerance
                warning('bounds:radius', 'The bounds for optimisation have restricted the optimal radius. Try increasing the bounds margin.');
            end
            if (azim_n < 2*step_tolerance && azim_LB > -pi) || (azim_n > 1 - 2*step_tolerance && azim_UB < pi) || ...
               (elev_n < 2*step_tolerance && elev_LB > 0) || (elev_n > 1 - 2*step_tolerance && elev_UB < pi/2)                                                      % The bounds are allowed to be hit if the angles are the spherical angle bounds
                warning('bounds:axis', 'The bounds for optimisation have restricted the optimal cylinder axis direction. Try increasing the bounds margin.');
            end

            % Check if the objective function is convex at the solution (i.e. if all eigenvalues of the Hessian are positive)
            % Note that as the optimising algorithm is quasi-Newton, convexity is not guaranteed
            Hessian_n_eigenvalues       = eig(Optimum_Hessian_matrix_n);
            Hessian_n_eigenvalues_sign  = sign(Hessian_n_eigenvalues);

            if sum(Hessian_n_eigenvalues_sign) < num_optim_geo_parameters
                warning('optimiser:convexity', 'The objective function is not convex at the found solution. The lowest Hessian eigenvalue is %.3g \n', min(Hessian_n_eigenvalues));
            end

            %--% Check for convergence %--%
            % If the difference between the current and previous geometry vectors is less than the threshold, it is said to have converged
            geometry_difference             = norm(inf_cylinder_geometry_prev_n - inf_cylinder_geometry_n);
            geometry_difference_list(iter)  = geometry_difference;

            if geometry_difference < UU_convergence_threshold
                Convergence = true;

            % Otherwise the geometry vector is updated
            else
                inf_cylinder_geometry_prev_n = inf_cylinder_geometry_n;
            end
        end

        % Final optimal geometry
        [Cylinder_radius, ~, Cylinder_direction, ~, Cylinder_centre]    = Geometry_Retrieval(inf_cylinder_geometry_n, geometry_LB, geometry_UB, gimbal_rotation_matrix, point_cloud_centroid, Scanner_loc_cell_r, Geometry_Indices);
        Optimal_Geometry                                                = struct('Cylinder_centre', Cylinder_centre, 'Cylinder_direction', Cylinder_direction, 'Cylinder_radius', Cylinder_radius);

        % Number of update iterations it took
        number_iterations = iter;

        % Removing excess entries 
        inf_cylinder_geometry_n_matrix(iter + 1: max_UU_iterations, :) = [];
        avg_Sigma_det_list(iter + 1 : max_UU_iterations)               = [];
        number_optim_steps_list(iter + 1: max_UU_iterations)           = [];
        geometry_difference_list(iter + 1: max_UU_iterations)          = [];

    %% Computation of the uncertainty %%
        % The range and radial uncertainty of each point
        [sigma_radial_cell, sigma_prop_cell, ~, ~] = Cylindrical_Object_Uncertainty(Cylinder_centre, Cylinder_direction, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);

        % The 3D distributions
        alpha                       = 1 - Confidence_interval/100;          % The confidence interval is changed to alpha
        Distr_Diagnostics           = false;
        Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

    %% Optimiser diagnostics information %%
        % Note that they are the final optimiser's diagnostics

        % The gradient and Hessian are unnormalised
        grad_delta_list         = geometry_UB - geometry_LB;
        Optimum_gradient_list   = grad_delta_list .* Optimum_gradient_list_n';

        hess_delta_matrix                   = grad_delta_list' * grad_delta_list;
        Optimum_Hessian_matrix              = hess_delta_matrix .* Optimum_Hessian_matrix_n;
        Optimum_second_derivative_list      = diag(Optimum_Hessian_matrix);
        Optimum_second_derivative_list_n    = diag(Optimum_Hessian_matrix_n);

        % For clarity later, the information is stored including the variable names for the first and second derivatives
        Optimum_First_Derivatives   = struct('first_derivative_list', Optimum_gradient_list);
        Norm_Opt_First_Derivatives  = struct('first_derivative_list', Optimum_gradient_list_n);
        Optimum_Second_Derivatives  = struct('second_derivative_list', Optimum_second_derivative_list);
        Norm_Opt_Second_Derivatives = struct('second_derivative_list', Optimum_second_derivative_list_n);
        Optimum_Hessian             = struct('Hessian_matrix', Optimum_Hessian_matrix);
        Norm_Optimum_Hessian        = struct('Hessian_matrix', Optimum_Hessian_matrix_n);

        geom_parameter_labels   = {'Circle_centre_a', 'Circle_centre_b', 'Azim_angle', 'Elev_angle', 'Radius'};
        geom_parameter_units    = {'m', 'm', 'rad', 'rad', 'm'};
        geom_parameter_indices  = [centre_ind, azim_ind, elev_ind, radius_ind];

        for i = 1 : num_optim_geo_parameters
            % The i'th parameter's label
            parameter_label_i   = geom_parameter_labels{i};
            parameter_ind_i     = geom_parameter_indices(i);

            % The first and second derivatives
            Optimum_First_Derivatives.(parameter_label_i)   = Optimum_gradient_list(parameter_ind_i);
            Norm_Opt_First_Derivatives.(parameter_label_i)  = Optimum_gradient_list_n(parameter_ind_i);
            Optimum_Second_Derivatives.(parameter_label_i)  = Optimum_second_derivative_list(parameter_ind_i);
            Norm_Opt_Second_Derivatives.(parameter_label_i) = Optimum_second_derivative_list_n(parameter_ind_i);

            for j = i : num_optim_geo_parameters
                % The j'th parameter's label
                parameter_label_j   = geom_parameter_labels{j};
                parameter_ind_j     = geom_parameter_indices(j);

                % The second derivative between the i'th and j'th parameter (which is saved twice for symmetry)
                Optimum_Hessian.(parameter_label_i).(parameter_label_j).cross_derivative        = Optimum_Hessian_matrix(parameter_ind_i, parameter_ind_j);
                Optimum_Hessian.(parameter_label_j).(parameter_label_i).cross_derivative        = Optimum_Hessian_matrix(parameter_ind_i, parameter_ind_j);
                Norm_Optimum_Hessian.(parameter_label_i).(parameter_label_j).cross_derivative   = Optimum_Hessian_matrix_n(parameter_ind_i, parameter_ind_j);
                Norm_Optimum_Hessian.(parameter_label_j).(parameter_label_i).cross_derivative   = Optimum_Hessian_matrix_n(parameter_ind_i, parameter_ind_j);
            end
        end
        
    %% Optimiser steps %%
        % The uncertainty updates are combined
        Infinite_cylinder_geometry_steps_n  = vertcat(Infinite_cylinder_geometry_steps_cell_n{:});
        Gradient_steps_n                    = vertcat(Gradient_steps_cell_n{:});
        Objective_value_steps_n             = vertcat(Objective_value_steps_cell_n{:});
    
        % Step sizes (based on the normalised geometry)
        Geometry_steps_n_delta_matrix       = diff(Infinite_cylinder_geometry_steps_n, 1, 1);
        Optimiser_step_size_list            = [0; sqrt(sum(Geometry_steps_n_delta_matrix.^2, 2))];              % A 0 value is appended for the first step

        % The optimisation steps are non-normalised
        Infinite_cylinder_geometry_steps    = Infinite_cylinder_geometry_steps_n .* (geometry_UB - geometry_LB) + geometry_LB;
        Gradient_steps                      = Gradient_steps_n .* grad_delta_list;
        Objective_value_steps               = Objective_value_steps_n;                                          % Matlab doesn't like global variables as outputs
        number_optimiser_steps_total        = sum(number_optim_steps_list);
        
        % A structure is created with this and other diagnostics data
        Geometry_Steps  = struct('Circle_centre_a', Infinite_cylinder_geometry_steps(:, centre_ind(1)), 'Circle_centre_b', Infinite_cylinder_geometry_steps(:, centre_ind(2)), 'Azim_angle', Infinite_cylinder_geometry_steps(:, azim_ind), 'Elev_angle', Infinite_cylinder_geometry_steps(:, elev_ind), 'Radius', Infinite_cylinder_geometry_steps(:, radius_ind));
        Gradient_Steps  = struct('Circle_centre_a', Gradient_steps(:, centre_ind(1)), 'Circle_centre_b', Gradient_steps(:, centre_ind(2)), 'Azim_angle', Gradient_steps(:, azim_ind), 'Elev_angle', Gradient_steps(:, elev_ind), 'Radius', Gradient_steps(:, radius_ind));
        Optimiser_Steps = struct('Geometry_Steps', Geometry_Steps, 'Gradient_Steps', Gradient_Steps, 'Step_sizes', Optimiser_step_size_list, 'Objective_value_steps', Objective_value_steps, 'num_optimiser_steps', number_optimiser_steps_total);

        Optimiser_Diagnostics = struct('parameter_labels', {geom_parameter_labels}, 'parameter_units', {geom_parameter_units}, 'num_parameters', num_optim_geo_parameters, 'Optimum_First_Derivatives', Optimum_First_Derivatives, 'Norm_Opt_First_Derivatives', Norm_Opt_First_Derivatives, 'Optimum_Second_Derivatives', Optimum_Second_Derivatives, 'Norm_Opt_Second_Derivatives', Norm_Opt_Second_Derivatives, 'Optimum_Hessian', Optimum_Hessian, 'Norm_Optimum_Hessian', Norm_Optimum_Hessian, 'Optimiser_Steps', Optimiser_Steps, 'number_UU_iterations', number_iterations);
        
    %% Fitting information %%
        % Printed statement
        if Print == true
            % The objective value is converted into the expected Mahalanobis distance
            Expected_Mahal_distance = Objective_value * Objective_value_init;

            fprintf('The infinite cylinder has been fitted with O =  %.4g, E[M^2] = %.4g \n', Objective_value, Expected_Mahal_distance);
            fprintf('   Radius:     %.3g m \n',             Cylinder_radius);
            fprintf('   Centre:     [%.3g, %.3g, %.3g] \n', Cylinder_centre);
            fprintf('   Direction:  [%.3g, %.3g, %.3g] \n', Cylinder_direction);
            fprintf('It took a total of %i steps over %i loops. \n', number_optimiser_steps_total, number_iterations);
        end      
        
        % Plots showing the point cloud, fitted cylinder and optimisation history
        if Plot == true
            %--% Geometry overview %--%
            % The number of coordinates in each ellipsoid and the cylinder
            number_coord = 1e3;     
            
            % Colours for the point cloud of each scanner
            scanner_colours = cbrewer('qual', 'Set2', max(number_scanners, 3));     % The colours used for the scanners
            scanner_colours = max(scanner_colours, 0);
            scanner_colours = min(scanner_colours, 1);
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
                
            % The initial cylinder surface
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius_init, cylinder_length, cylinder_centre_init, cylinder_direction_init, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial infinite cylinder');
            
            % The optimised cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(Cylinder_radius, cylinder_length, Cylinder_centre, Cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Optimised infinite cylinder');

            % The point cloud and its uncertainty
            distribution_mu_cell    = Point_Cloud_Distributions.distribution_mu_cell;
            distribution_radii_cell = Point_Cloud_Distributions.distribution_radii_cell;
            distribution_axes_cell  = Point_Cloud_Distributions.distribution_axes_cell;
            
            number_distributions_list   = Point_Cloud_Distributions.number_distributions_list;
            num_cumulative_distr_list   = [0, cumsum(number_distributions_list)];
            
            for s = 1 : number_scanners
                %--% Scanner vectors %--%
                % Distance from scanner to cylinder
                scanner_loc     = Scanner_loc_cell{s};
                scanner_vector  = scanner_loc - Cylinder_centre;
                range           = sqrt(sum(scanner_vector.^2));
                
                % The vector
                scanner_vector_scaled   = cylinder_length * scanner_vector / norm(scanner_vector);
                scanner_string          = sprintf('Scanner %g, R = %.3g m', s, range);
                scanner_colour          = scanner_colours(s, :);
                plot3(Cylinder_centre(1) + [0, scanner_vector_scaled(1)], Cylinder_centre(2) + [0, scanner_vector_scaled(2)], Cylinder_centre(3) + [0, scanner_vector_scaled(3)], 'LineWidth', 2, 'color', scanner_colour, 'DisplayName', scanner_string);
                
                %--% Scanner's point cloud %--%
                start_ind   = num_cumulative_distr_list(s) + 1;
                end_ind     = num_cumulative_distr_list(s + 1);
                
                for d = start_ind : end_ind
                    % Distribution properties
                    distribution_mu     = distribution_mu_cell{d};
                    distribution_radii  = distribution_radii_cell{d};
                    distribution_axes   = distribution_axes_cell{d};

                    % Surface
                    [ellipsoid_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_coord);

                    x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                    y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                    z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                    surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', scanner_colour, 'FaceAlpha', 0.10, 'DisplayName', sprintf('Point cloud %g uncertainty, \\alpha = %.3g', s, alpha));

                    % Mu
                    sc_mu = scatter3(distribution_mu(1), distribution_mu(2), distribution_mu(3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', scanner_colour, 'DisplayName', '\mu');

                    if d > start_ind
                        sc_mu.HandleVisibility      = 'Off';
                        surf_el.HandleVisibility    = 'Off';
                    end
                end 
            end
            
            % Aspect ratio
            axis equal

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')
            
            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off
                        
            %--% Optimisation parameter history %--%
            optim_parameter_labels  = {'Q', '\Delta'};
            optim_parameter_units   = {'-', '-'};
            num_optim_parameters    = length(optim_parameter_labels);
            
            optim_parameter_colours = cbrewer('qual', 'Set2', max(num_optim_parameters, 3));
            optim_parameter_colours = max(optim_parameter_colours, 0);
            optim_parameter_colours = min(optim_parameter_colours, 1);
            
            optim_parameter_history_matrix = [Objective_value_steps, Optimiser_step_size_list];
            
            figure(2)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   
            
            for q = 1 : num_optim_parameters
                % This parameter's data
                optim_parameter_label   = optim_parameter_labels{q};
                optim_parameter_unit    = optim_parameter_units{q};
                optim_parameter_colour  = optim_parameter_colours(q, :);
                optim_parameter_history = optim_parameter_history_matrix(:, q);
                
                % The subplot
                subplot(1, num_optim_parameters, q)
                hold on
                grid on

                plot(1 : number_optimiser_steps, optim_parameter_history, 'LineWidth', 2, 'Color', optim_parameter_colour);

                % Axes
                xlim([1, number_optimiser_steps]);
                xlabel('Optim. step [-]');
                ylabel(sprintf('%s [%s]', optim_parameter_label, optim_parameter_unit));

                % Font
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off   
            end
            
            %--% Geometry parameter (gradient) history %--%    
            % Colours
            geometry_colours    = cbrewer('qual', 'Set1', num_optim_geo_parameters);
            geometry_colours    = max(geometry_colours, 0);
            geometry_colours    = min(geometry_colours, 1);

            for g = 1 : num_optim_geo_parameters
                % This geometry type's information
                geometry_colour     = geometry_colours(g, :);
                geometry_label      = geom_parameter_labels{g};
                geometry_unit       = geom_parameter_units{g};

                geometry_string     = strrep(geometry_label, '_', ' ');

                geometry_steps          = Geometry_Steps.(geometry_label);
                geometry_gradient_steps = Gradient_Steps.(geometry_label);

                % The plot
                figure(10 + g)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])           

                % Geometry variable
                subplot(1, 2, 1)
                hold on
                grid on

                plot(1 : number_optimiser_steps, geometry_steps, 'LineWidth', 2, 'color', geometry_colour);

                % Axes
                xlim([1, number_optimiser_steps]);
                xlabel('Optim. step [-]');
                ylabel(sprintf('%s [%s]', geometry_string, geometry_unit));

                % Font
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                % Its gradient
                subplot(1, 2, 2)
                hold on
                grid on

                plot(1 : number_optimiser_steps, geometry_gradient_steps, 'LineWidth', 2, 'color', geometry_colour);

                % Axes
                xlim([1, number_optimiser_steps]);
                xlabel('Optim. step [-]');
                ylabel(sprintf('dQ/d%s [1/%s]', geometry_string, geometry_unit));

                % Font
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off   
            end

            %--% Iterative shape fitting diagnostics %--%
            % Geometry parameter updates
            inf_cylinder_geometry_matrix = inf_cylinder_geometry_n_matrix .* (geometry_UB - geometry_LB) + geometry_LB;

            figure(100)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   

            for g = 1 : num_optim_geo_parameters
                % This geometry type's information
                geometry_list       = inf_cylinder_geometry_matrix(:, g);
                geometry_colour     = geometry_colours(g, :);
                geometry_label      = geom_parameter_labels{g};
                geometry_unit       = geom_parameter_units{g};
    
                geometry_string     = strrep(geometry_label, '_', ' ');
 
                % The subplot
                subplot(2, ceil(num_optim_geo_parameters / 2), g)
                hold on
                grid on

                plot(geometry_list, 'LineWidth', 2, 'color', geometry_colour)

                xlabel('Iteration [-]');
                ylabel(sprintf('%s [%s]', geometry_string, geometry_unit))

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off   
            end

            % Diagnostics
            iter_diagnostics_cell   = {avg_Sigma_det_list, number_optim_steps_list, geometry_difference_list};
            iter_diag_label_cell    = {sprintf('avg. |%s|', '\Sigma'), 'number steps', sprintf('%s g', '\Delta')};
            iter_diag_unit_cell     = {'m^6', '-', '-'};
            number_iter_diag_sets   = length(iter_diag_unit_cell);

            iter_diag_colours = cbrewer('qual', 'Set1', max(3, number_iter_diag_sets));
            iter_diag_colours = max(iter_diag_colours, 0);
            iter_diag_colours = min(iter_diag_colours, 1);

            figure(101)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   

            for i = 1 : number_iter_diag_sets
                % This diagnostic's sets info
                iter_diagnostics_list   = iter_diagnostics_cell{i};
                iter_diag_label         = iter_diag_label_cell{i};
                iter_diag_unit          = iter_diag_unit_cell{i};
                iter_diag_colour        = iter_diag_colours(i, :);

                subplot(1, number_iter_diag_sets, i);
                hold on
                grid on

                plot(iter_diagnostics_list, 'LineWidth', 2, 'color', iter_diag_colour);

                xlabel('Iteration [-]');
                ylabel(sprintf('%s [%s]', iter_diag_label, iter_diag_unit));
                ylim([0, 1.2*max(iter_diagnostics_list)]);

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off   
            end

            %--% Pause message %--%
            disp('The script will continue and the plots will close when a key is pressed');
            pause();

            close all;             
        end

    %% Local functions %%
        % Geometry variable retrieval from the vector
        function [cylinder_radius, cylinder_direction_r, cylinder_direction, cylinder_centre_r, cylinder_centre] = Geometry_Retrieval(inf_cylinder_geometry_n, geometry_LB, geometry_UB, gimbal_rotation_matrix, point_cloud_centroid, Scanner_loc_cell_r, Geometry_Indices)
            % Un-normalisation
            inf_cylinder_geometry = inf_cylinder_geometry_n .* (geometry_UB - geometry_LB) + geometry_LB;
            
            % Geometry indices
            [centre_ind, azim_ind, elev_ind, radius_ind] = deal(Geometry_Indices.centre, Geometry_Indices.azim_angle, Geometry_Indices.elev_angle, Geometry_Indices.radius);

            % The spherical angles and resulting cylinder axis
            [azim_angle_r, elev_angle_r]    = deal(inf_cylinder_geometry(azim_ind), inf_cylinder_geometry(elev_ind));
            [cylinder_direction_r, ~, ~]    = Vector_Spherical_Angle_Conversion([], azim_angle_r, elev_angle_r);
    
            cylinder_direction              = (gimbal_rotation_matrix' * cylinder_direction_r')';
    
            % The circle centre is converted into the cylinder centre
            circle_centre_r                 = inf_cylinder_geometry(centre_ind);
            height_r                        = 0;
            [~, ~, cylinder_centre_r, ~]    = Circle_Cylinder_Centre_Conversion(circle_centre_r, height_r, [], cylinder_direction_r, origin, Scanner_loc_cell_r);
    
            cylinder_centre_c               = (gimbal_rotation_matrix' * cylinder_centre_r')';
            cylinder_centre                 = cylinder_centre_c + point_cloud_centroid;
    
            % The cylinder radius equals the circular cross-section's radius
            cylinder_radius                 = inf_cylinder_geometry(Geometry_Indices.radius);
        end

        % Optimisation path function
        function stop = Optim_Path(inf_cylinder_geometry_n, Optimisation_Values, ~)
            stop = false;

            % Concatenate objective value, geometry and various other terms computed during each iteration
            Infinite_cylinder_geometry_steps_cell_n{iter}   = [Infinite_cylinder_geometry_steps_cell_n{iter}; inf_cylinder_geometry_n];
            Gradient_steps_cell_n{iter}                     = [Gradient_steps_cell_n{iter}; Optimisation_Values.gradient'];
            Objective_value_steps_cell_n{iter}              = [Objective_value_steps_cell_n{iter}; Optimisation_Values.fval];
        end
end