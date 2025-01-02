% The optimal cylindrical cross-section and axis are determined with the analytical Hessian and gradient (for which constant uncertainty is assumed)
% As objective E[M^2] is used with the line approximation. The history of the step size and geometry are also given as outputs
% To improve behaviour of the algorithm the search space can be bounded and for the update step the Hessian and gradient normalised using these bounds

% Note: Inputs are expected to be centered at the point cloud centroid (i.e. it is taken to be the origin in this script)

function [Optimal_Inf_Cyl_Geometry_c, Optimiser_Diagnostics, Objective_value] = Iterative_Infinite_Cylinder_Fitting(Init_Cyl_Geometry_c, Fitting_Parameters, Point_Cloud_Distributions_c, Scanning_Parameters_c, Statistical_Values, Output_Decisions)
    
    %% Structure inputs %%
        % Cylinder geometry
        init_cylinder_centre_c  = Init_Cyl_Geometry_c.Cylinder_centre;
        init_cylinder_dir       = Init_Cyl_Geometry_c.Cylinder_direction;
        init_cylinder_radius    = Init_Cyl_Geometry_c.Cylinder_radius;
        cylinder_length         = Init_Cyl_Geometry_c.Cylinder_length;
       
        % Fitting parameters
        distance_moment         = Fitting_Parameters.distance_moment;
        Point_weights_list      = Fitting_Parameters.point_weights_list;
        bounds_margin           = Fitting_Parameters.bounds_margin;
        Distance_Computation    = Fitting_Parameters.Distance_Computation;

        % Scanner locations
        number_scanners         = Scanning_Parameters_c.number_scanners;
        Scanner_loc_cell_c      = Scanning_Parameters_c.Scanner_loc_cell;

        % Output decisions
        Print                   = Output_Decisions.Print;
        Plot                    = Output_Decisions.Plot;
  
    %% Manual inputs %%
        % Central difference step size factor for the gradient and Hessian
        step_size_factor        = 1e-4;             % [-] Multiplied by the difference in geometry variable bounds

        % Diagnostic outputs
        Mahal_Diagnostics       = false;            % [true, false]
        Minimiser_Print         = true;             % [true, false]
        Minimiser_Plot          = false;            % [true, false]

    %% Geometry indices %%
        % To ensure consistency of the location for geometrical variables in an array
        num_dim = length(init_cylinder_dir);
        
        centre_ind      = 1 : num_dim - 1;
        radius_ind      = max(centre_ind) + 1;
        azim_ind        = radius_ind + 1;
        elev_ind        = azim_ind + 1;
        num_geom_var    = max([centre_ind, radius_ind, azim_ind, elev_ind]);

    %% Solution bounds %%
        % Reasonable bounds for the optimisation process are based on the initial geometry estimate for the circle centre and radius
        origin                          = zeros(1, num_dim);
        [init_circle_centre_c, ~, ~, ~] = Circle_Cylinder_Centre_Conversion([], [], init_cylinder_centre_c, init_cylinder_dir, origin, Scanner_loc_cell_c);

        circle_centre_LB    = init_circle_centre_c - bounds_margin * init_cylinder_radius;
        circle_centre_UB    = init_circle_centre_c + bounds_margin * init_cylinder_radius;
        
        radius_LB   = init_cylinder_radius / (1 + bounds_margin);
        radius_UB   = init_cylinder_radius * (1 + bounds_margin);

        % The direction vectors can deviate by the radius times the margin, at the cylinder length's distance
        % Note that this deviation only results in allowed changes to the elevation angle. The azimuth angle remains unbounded
        [~, ~, ~, ~, azim_angle, elev_angle] = Vector_Spherical_Angle_Conversion(init_cylinder_dir, [], []);

        elev_deviation      = atan(init_cylinder_radius * bounds_margin / cylinder_length);
        elev_LB             = max(-pi/2, elev_angle - elev_deviation);
        elev_UB             = min(pi, elev_angle + elev_deviation);
        [azim_LB, azim_UB]  = deal(-pi, pi);

        % The combined sets of bounds
        geometry_LB                                                 = zeros(1, num_geom_var);
        geometry_LB([centre_ind, azim_ind, elev_ind, radius_ind])   = [circle_centre_LB, azim_LB, elev_LB, radius_LB];

        geometry_UB                                                 = ones(1, num_geom_var);
        geometry_UB([centre_ind, azim_ind, elev_ind, radius_ind])   = [circle_centre_UB, azim_UB, elev_UB, radius_UB];

        % Resulting step sizes
        step_size_list  = step_size_factor * (geometry_UB - geometry_LB);

    %% Shape fitting %%
        % Initial geometry estimate
        init_inf_cyl_geometry                                               = zeros(1, num_geom_var);
        init_inf_cyl_geometry([centre_ind, radius_ind, azim_ind, elev_ind]) = [init_circle_centre_c, init_cylinder_radius, azim_angle, elev_angle];

        % Its objective value for normalisation
        Objective_value_init = Objective_Function(init_inf_cyl_geometry, 1, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);
           
        % Function handles for the objective, gradient and Hessian
        Objective_fun   = @(inf_cyl_geometry) Objective_Function(inf_cyl_geometry, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);
        Gradient_fun    = @(inf_cyl_geometry) Gradient_Function(inf_cyl_geometry, Objective_value_init, step_size_list, num_geom_var, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);
        Hessian_fun     = @(inf_cyl_geometry) Hessian_Function(inf_cyl_geometry, Objective_value_init, step_size_list, num_geom_var, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);
            
        % Steepest descent minimisation which can continue after hitting the bounds and from hitting a saddle point
        [opt_inf_cyl_geometry, Objective_value, Optimiser_Steps] = Steepest_Descent_Minimiser(init_inf_cyl_geometry, geometry_LB, geometry_UB, Objective_fun, Gradient_fun, Hessian_fun, Minimiser_Print, Minimiser_Plot);

        % The geometry vector is turned into geometry parameters
        [Optimal_Inf_Cyl_Geometry_c, ~, Cylinder_centre_c, Cylinder_radius, Cylinder_dir] = Geometry_Retrieval(opt_inf_cyl_geometry);

        % Gradient and Hessian values at the optimum
        Opt_gradient_list   = Gradient_fun(opt_inf_cyl_geometry);
        Opt_Hessian_matrix  = Hessian_fun(opt_inf_cyl_geometry);
        
    %% Optimiser steps %%
        %--% Result structures %--%
        % Data for specific variables is added and a generic diagnostics structure is created
        geom_variable_labels    = {'Circle_centre_a', 'Circle_centre_b', 'Azim_angle', 'Elev_angle', 'Radius'};
        geom_variable_units     = {'m', 'm', 'rad', 'rad', 'm'};
        geom_variable_indices   = [centre_ind, azim_ind, elev_ind, radius_ind];

        Optimiser_Diagnostics   = struct('parameter_labels', {geom_variable_labels}, 'parameter_units', {geom_variable_units}, 'num_parameters', num_geom_var, ...
                                         'Optimum_First_Derivatives', struct('first_derivative_list', Opt_gradient_list), 'Optimum_Second_Derivatives', struct('second_derivative_list', diag(Opt_Hessian_matrix)'), 'Optimum_Hessian', struct('Hessian_matrix', Opt_Hessian_matrix), ...
                                         'Norm_Opt_First_Derivatives', struct('first_derivative_list', NaN(1, num_geom_var)), 'Norm_Opt_Second_Derivatives', struct('second_derivative_list', NaN(1, num_geom_var)), 'Norm_Optimum_Hessian', struct('Hessian_matrix', NaN(num_geom_var)));          % Note that no normalisation took place

        num_optimiser_steps     = Optimiser_Steps.num_optimiser_steps;

        for i = 1 : num_geom_var
            % This variable's name and index
            variable_label_i    = geom_variable_labels{i};
            variable_ind_i      = geom_variable_indices(i);

            % Values at each step
            Optimiser_Steps.Geometry_Steps.(variable_label_i)           = Optimiser_Steps.Variable_steps(:, variable_ind_i);
            Optimiser_Steps.Gradient_Steps.(variable_label_i)           = Optimiser_Steps.Gradient_steps(:, variable_ind_i);
            Optimiser_Steps.Second_Derivative_Steps.(variable_label_i)  = Optimiser_Steps.Second_derivative_steps(:, variable_ind_i);                          

            % Values at the optimum
            Optimiser_Diagnostics.Optimum_First_Derivatives.(variable_label_i)      = Optimiser_Steps.Gradient_steps(num_optimiser_steps, variable_ind_i);
            Optimiser_Diagnostics.Norm_Opt_First_Derivatives.(variable_label_i)     = NaN(num_optimiser_steps, 1);                                                      % Normalisation was not applied
            Optimiser_Diagnostics.Optimum_Second_Derivatives.(variable_label_i)     = Optimiser_Steps.Second_derivative_steps(num_optimiser_steps, variable_ind_i);
            Optimiser_Diagnostics.Norm_Opt_Second_Derivatives.(variable_label_i)    = NaN(num_optimiser_steps, 1);                                                      % Normalisation was not applied

            % Hessian cross-derivatives
            for j = i : num_geom_var
                % The j'th variable's name and index
                variable_label_j    = geom_variable_labels{j};
                variable_ind_j      = geom_variable_indices(j);

                % The second derivative between the i'th and j'th variables (which is saved twice for symmetry)
                Optimiser_Diagnostics.Optimum_Hessian.(variable_label_i).(variable_label_j).cross_derivative        = Optimiser_Steps.Hessian_steps(variable_ind_i, variable_ind_j, num_optimiser_steps);
                Optimiser_Diagnostics.Optimum_Hessian.(variable_label_j).(variable_label_i).cross_derivative        = Optimiser_Steps.Hessian_steps(variable_ind_j, variable_ind_i, num_optimiser_steps);
                Optimiser_Diagnostics.Norm_Optimum_Hessian.(variable_label_j).(variable_label_i).cross_derivative   = NaN;                                              % Normalisation was not applied
                Optimiser_Diagnostics.Norm_Optimum_Hessian.(variable_label_i).(variable_label_j).cross_derivative   = NaN; 
            end
        end

        % The optimiser steps are added to the diagnostics structure, without superfluous fields
        Optimiser_Steps                         = rmfield(Optimiser_Steps, {'Variable_steps', 'Gradient_steps', 'Search_direction_steps', 'Second_derivative_steps'});
        Optimiser_Diagnostics.Optimiser_Steps   = Optimiser_Steps;

    %% Printed messages %%
        % Final E[M^2] and geometry
        if Print == true
            fprintf('The infinite cylinder has been fitted in %i steps with E[M^2] = %.4g \n', num_optimiser_steps, Objective_value);
            fprintf('   Radius:     %.3g m \n',             Cylinder_radius);
            fprintf('   Centre:     [%.3g, %.3g, %.3g] \n', Cylinder_centre_c);
            fprintf('   Direction:  [%.3g, %.3g, %.3g] \n', Cylinder_dir);
        end    

    %% Plots %%
        % Plot showing the point cloud with fitted cylinder
        if Plot == true
            % The number of coordinates in each ellipsoid and the cylinder
            number_coord = 1e3;     

            % The colours used for the scanners
            scanner_colours = cbrewer('qual', 'Set2', max(number_scanners, 3));     
            scanner_colours = max(scanner_colours, 0);
            scanner_colours = min(scanner_colours, 1);
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
                
            % The initial cylinder surface
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~] = Cylinder_Surface_Generator(init_cylinder_radius, cylinder_length, init_cylinder_centre_c, init_cylinder_dir, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'm', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial infinite cylinder');
            
            % The optimised cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(Cylinder_radius, cylinder_length, Cylinder_centre_c, Cylinder_dir, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Optimised infinite cylinder');

            % The point cloud and its uncertainty
            distribution_mu_cell    = Point_Cloud_Distributions_c.distribution_mu_cell;
            distribution_radii_cell = Point_Cloud_Distributions_c.distribution_radii_cell;
            distribution_axes_cell  = Point_Cloud_Distributions_c.distribution_axes_cell;
            distr_alpha             = Point_Cloud_Distributions_c.alpha;

            number_distributions_list   = Point_Cloud_Distributions_c.number_distributions_list;
            num_cumulative_distr_list   = [0, cumsum(number_distributions_list)];
            
            for s = 1 : number_scanners
                %--% Scanner vectors
                % Distance from scanner to cylinder
                scanner_loc     = Scanner_loc_cell_c{s};
                scanner_vector  = scanner_loc - Cylinder_centre_c;
                range           = sqrt(sum(scanner_vector.^2));
                
                % The vector
                scanner_vector_scaled   = cylinder_length * scanner_vector / norm(scanner_vector);
                scanner_string          = sprintf('Scanner %g, R = %.3g m', s, range);
                scanner_colour          = scanner_colours(s, :);
                plot3(Cylinder_centre_c(1) + [0, scanner_vector_scaled(1)], Cylinder_centre_c(2) + [0, scanner_vector_scaled(2)], Cylinder_centre_c(3) + [0, scanner_vector_scaled(3)], 'LineWidth', 2, 'color', scanner_colour, 'DisplayName', scanner_string);
                
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

                    surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', scanner_colour, 'FaceAlpha', 0.10, 'DisplayName', sprintf('Point cloud %g uncertainty, \\alpha = %.3g', s, distr_alpha));

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

            % Pause message
            disp('The optimal infinite cylinder has been fitted. The plot will close when a key is pressed');
            pause();

            close(1);             
        end
        
    %% Local functions %%
        % Conversion from geometry vector array to more descriptive geometry 
        function [Cylinder_Geometry, circle_centre_c, cylinder_centre_c, cylinder_radius, cylinder_dir] = Geometry_Retrieval(inf_cyl_geometry)
            % Vector to cylinder geometry
            [circle_centre_c, cylinder_radius, gamma, epsilon] = deal(inf_cyl_geometry(centre_ind), inf_cyl_geometry(radius_ind), inf_cyl_geometry(azim_ind), inf_cyl_geometry(elev_ind));
            
            [~, ~, ~, cylinder_dir, ~, ~]   = Vector_Spherical_Angle_Conversion([], gamma, epsilon);
            circle_height                   = 0;                                                                % Note that the height above the cross-section is irrelevant when the problem is projected onto it
            [~, ~, ~, cylinder_centre_c]    = Circle_Cylinder_Centre_Conversion(circle_centre_c, circle_height, [], cylinder_dir, origin, Scanner_loc_cell_c);

            % Structure
            Cylinder_Geometry = struct('Cylinder_centre', cylinder_centre_c, 'Cylinder_radius', cylinder_radius, 'Cylinder_direction', cylinder_dir, 'Cylinder_length', cylinder_length);
        end

        % Objective function
        function Objective_value = Objective_Function(inf_cyl_geometry, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics)
            % The geometry
            [~, ~, cyl_centre_c, cyl_radius, cyl_dir] = Geometry_Retrieval(inf_cyl_geometry);

            % Expected squared Mahalanobis distance
            if strcmp(Distance_Computation, 'Line_Approx')
                [~, expected_Mahal_distance_list] = Expected_Mahal_Distance_Cylinder_Line_Approx(cyl_radius, cyl_centre_c, cyl_dir, Point_Cloud_Distributions_c, Scanner_loc_cell_c, distance_moment, Mahal_Diagnostics);
            elseif strcmp(Distance_Computation, 'Numerical')
                [~, expected_Mahal_distance_list] = Expected_Mahalanobis_Distance_Cylinder_Numerical(cyl_centre_c, cyl_radius, cyl_dir, distance_moment, Point_Cloud_Distributions_c, Statistical_Values, Mahal_Diagnostics);
            end
            
            % The weighted average is used as the objective
            expected_Mahal_distance_list_w  = Point_weights_list .* expected_Mahal_distance_list;
            Objective_value                 = mean(expected_Mahal_distance_list_w) / Objective_value_init;
        end

        % Gradient function 
        function gradient_vector = Gradient_Function(inf_cyl_geometry, Objective_value_init, step_size_list, num_geometry_variables, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics)
            % Gradient for each geometry variable
            gradient_vector = zeros(1, num_geometry_variables);

            for v = 1 : num_geometry_variables
                % This variable's step size
                step_size = step_size_list(v);
            
                % The objective value with a forward step for this dimension
                inf_cyl_geometry_plus       = inf_cyl_geometry;
                inf_cyl_geometry_plus(v)    = inf_cyl_geometry_plus(v) + step_size;
                objective_value_plus        = Objective_Function(inf_cyl_geometry_plus, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);
                
                % The objective value with a backward step for this dimension
                inf_cyl_geometry_min        = inf_cyl_geometry;
                inf_cyl_geometry_min(v)     = inf_cyl_geometry_min(v) - step_size;
                objective_value_min         = Objective_Function(inf_cyl_geometry_min, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);
                
                % The resulting gradient
                gradient            = (objective_value_plus - objective_value_min) / (2*step_size);
                gradient_vector(v)  = gradient;
            end
        end

        % Hessian function 
        function Hessian_matrix = Hessian_Function(inf_cyl_geometry, Objective_value_init, step_size_list, num_geometry_variables, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics)
            % The objective value at the current geometry estimate
            objective_value = Objective_Function(inf_cyl_geometry, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);

            % Hessian matrix
            Hessian_matrix  = zeros(num_geometry_variables);

            for a = 1 : num_geometry_variables
                % This variable's step size
                step_size_a = step_size_list(a);
            
                % The objective value with a forward step for this dimension
                inf_cyl_geometry_a      = inf_cyl_geometry;
                inf_cyl_geometry_a(a)   = inf_cyl_geometry_a(a) + step_size_a;
                objective_value_a       = Objective_Function(inf_cyl_geometry_a, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);

                for b = a : num_geometry_variables
                    % This variable's step size
                    step_size_b = step_size_list(b);
                
                    % The objective value with a forward step for this dimension
                    inf_cyl_geometry_b      = inf_cyl_geometry;
                    inf_cyl_geometry_b(b)   = inf_cyl_geometry_b(b) + step_size_b;
                    objective_value_b       = Objective_Function(inf_cyl_geometry_b, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);

                    % The objective value with a forward step in both dimensions
                    inf_cyl_geometry_ab     = inf_cyl_geometry_a;
                    inf_cyl_geometry_ab(b)  = inf_cyl_geometry_ab(b) + step_size_b;
                    objective_value_ab      = Objective_Function(inf_cyl_geometry_ab, Objective_value_init, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Statistical_Values, Point_weights_list, Mahal_Diagnostics);

                    % Hessian value for variables a and b
                    Hessian_ab              = (objective_value_ab - objective_value_a - objective_value_b + objective_value) / (step_size_a * step_size_b);
                    Hessian_matrix(a, b)    = Hessian_ab;
                    Hessian_matrix(b, a)    = Hessian_ab;
                end
            end
        end
end