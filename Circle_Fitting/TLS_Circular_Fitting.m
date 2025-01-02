% This script attempts to determine the shape parameters (x, y, r) of a tree stem slice (2D), given TLS data
% Circle fitting can be done through (uncertainty weighted) algebraic and geometric least-squares, Pratt least-squares or optimisation with likelihood and/or expected distance

% Note that geometry can be given to filter for outliers and duplicate the points around it, if desired
% If Input_Geometry is empty, filtering will be performed if desired, however duplication will not
% Note that Input_Geometry is expected to be [x, y, r]

function [Tree_centre_fit, Tree_radius_fit, Objective_value, sigma_radial_cell, sigma_range_cell, laser_radius_cell, Tree_stem_geometry_steps, Objective_value_steps] = TLS_Circular_Fitting(x_points_cell, y_points_cell, Input_Geometry, Scanning_Parameters, Scanner_Parameters, Fitting_Parameters, Output_Decisions)
 
    %% Manual Inputs %%
        % As many fits might be made, printing can be set separately here
        Print                   = false;    % [true, false]
    
        % Optimisation depends on the initial estimate, which can be updated iteratively
        max_iterations          = 020;      % The maximum number of iterations until convergence
        convergence_threshold   = 005;      % The maximum percentage of change between the current iteration and the average of the previous iterations 

    %% Data inputs %%
    
        % Scanning parameters
        Scanner_loc_cell    = Scanning_Parameters.Scanner_loc_cell;
        number_scanners     = Scanning_Parameters.number_scanners;
        scanning_theta_cell = Scanning_Parameters.scanning_theta_cell;
        
        % Scanner parameters        
        laser_resolution    = Scanner_Parameters.laser_resolution * 1e-3;       % Rad
        beam_divergence     = Scanner_Parameters.beam_divergence * 1e-3;        % Rad
        max_incidence_angle = deg2rad(Scanner_Parameters.max_incidence_angle);  % Rad
        beam_exit_diameter  = Scanner_Parameters.beam_exit_diameter * 1e-3;     % m
        range_bias          = Scanner_Parameters.range_bias * 1e-3;             % m
        sigma_range_device  = Scanner_Parameters.sigma_range_device * 1e-3;     % m
        
        % Fitting parameters
        Fitting_method      = Fitting_Parameters.Fitting_method;
        Weighting           = Fitting_Parameters.Weighting;
        distance_moment     = Fitting_Parameters.distance_moment;
        objective_balance   = Fitting_Parameters.objective_balance;
        Optim_Convergence   = Fitting_Parameters.Optim_Convergence;
        Point_Duplication   = Fitting_Parameters.Point_Duplication;
        Point_Mirroring     = Fitting_Parameters.Point_Mirroring;
        bounds_margin       = Fitting_Parameters.bounds_margin;
        Filtering           = Fitting_Parameters.Filtering;
        Weighted_Filtering  = Fitting_Parameters.Weighted_Filtering;
        L_N_threshold       = Fitting_Parameters.L_N_threshold;
        
        % Output parameters
        Diagnostics         = Output_Decisions.Diagnostics;

    %% Point cloud pre-processing %%
        % Check if geometry was given or not
        if isempty(Input_Geometry)
            Geometry_Given      = false;
        else
            Geometry_Given      = true;
        end
    
        % Point cloud filtering
        if Filtering == true
            Filter_Plot = false;
            
            Weighted_Filtering_Init  = false;
            [x_f_points_cell, y_f_points_cell, ~, ~] = Kernel_Density_Estimation_Filtering(x_points_cell, y_points_cell, L_N_threshold, Weighted_Filtering_Init, [], [], Print, Filter_Plot);    
        else
            x_f_points_cell = x_points_cell;
            y_f_points_cell = y_points_cell;
        end
        
        % Point cloud duplication
        if Point_Duplication == true && Geometry_Given == true
            [x_fd_points_cell, y_fd_points_cell] = Point_Cloud_Duplication(Input_Geometry, x_f_points_cell, y_f_points_cell, number_scanners);   
        else
            x_fd_points_cell = x_f_points_cell;
            y_fd_points_cell = y_f_points_cell;
        end
                            
    %% The shape fitting constraints %%
        Constraint_Plot = false;
        [central_theta_list, delta_inner_list, delta_outer_list] = Bounds_from_Beams_2D(Scanner_loc_cell, scanning_theta_cell, x_points_cell, y_points_cell, laser_resolution, beam_divergence, Constraint_Plot);
  
    %% Shape fitting the point cloud %%
        %--% The best shape parameters are determined through optimisation %--%
        if strcmp(Fitting_method, 'Optim')          
            %--% Initial estimate of the tree stem parameters %--%
            if Geometry_Given == false
                Initial_Fitting_Method  = 'LS';
                Initial_Weighting       = true;
                [Tree_stem_x_init, Tree_stem_y_init, Tree_radius_init] = Least_Squares_Circle_Fitting(x_fd_points_cell, y_fd_points_cell, Initial_Fitting_Method, Initial_Weighting, Scanner_loc_cell, beam_divergence, sigma_range_device, range_bias, max_incidence_angle, Point_Duplication);
                
                Tree_stem_geometry_init = [Tree_stem_x_init, Tree_stem_y_init, Tree_radius_init];                        
            else
                Tree_stem_geometry_init = Input_Geometry;
            end

            % The objective value of the initial geometry
            [sigma_radial_cell, sigma_range_cell, ~, ~] = Circular_Object_Uncertainty(Tree_stem_geometry_init(1), Tree_stem_geometry_init(2), Tree_stem_geometry_init(3), Scanner_loc_cell, x_fd_points_cell, y_fd_points_cell, beam_divergence, sigma_range_device, max_incidence_angle, Point_Duplication);                        
            laser_radius_cell = Laser_Beamwidth(beam_exit_diameter, beam_divergence, x_fd_points_cell, y_fd_points_cell, Scanner_loc_cell, Point_Duplication);

            Diagnostics_Likelihood = false;
            [~, ~, Avg_max_likelihood, ~, ~, ~] = Circle_Fit_Likelihood(x_fd_points_cell, y_fd_points_cell, sigma_radial_cell, sigma_range_cell, laser_radius_cell, Tree_stem_geometry_init(1), Tree_stem_geometry_init(2), Tree_stem_geometry_init(3), Scanner_loc_cell, range_bias, Point_Duplication, Point_Mirroring, Diagnostics_Likelihood);

            Norm_expected_distance = 1;     % Note that it is 1, as there is no initial estimate to normalise it with
            
            [Geometry_Function_Handle, Geometry_variable_bounds] = Shape_Basis_Function();
            Objective_value_old = Objective_Function(Avg_max_likelihood, Norm_expected_distance, objective_balance, Geometry_Function_Handle, Tree_stem_geometry_init, Geometry_variable_bounds, central_theta_list, delta_inner_list, delta_outer_list, Scanner_loc_cell);

            %--% Iterative optimisation loop %--%
            if Optim_Convergence == false   % If convergence is not required, only one iteration is performed
                max_iterations = 1;
            end
            
            iter                = 0;
            convergence         = false;
            
            while convergence == false && iter <= max_iterations
                iter = iter + 1;
                 
                %--% Point cloud preprocessing %--%
                % Point cloud filtering
                if Filtering == true
                    [x_f_points_cell, y_f_points_cell, ~, ~] = Kernel_Density_Estimation_Filtering(x_points_cell, y_points_cell, L_N_threshold, Weighted_Filtering, sigma_range_cell, sigma_radial_cell, Print, Filter_Plot);
                else
                    x_f_points_cell = x_points_cell;
                    y_f_points_cell = y_points_cell;
                end

                % Point cloud duplication
                if Point_Duplication == true
                    [x_fd_points_cell, y_fd_points_cell] = Point_Cloud_Duplication(Tree_stem_geometry_init, x_f_points_cell, y_f_points_cell, number_scanners);   
                else
                    x_fd_points_cell = x_f_points_cell;
                    y_fd_points_cell = y_f_points_cell;
                end

                %--% Point cloud properties %--%
                % The uncertainty corresponding to this initial estimate is used to compute the likelihood/weight the samples (such that the optimiser is not tempted to change geometry to minimise uncertainty)
                [sigma_radial_cell, sigma_range_cell, ~, ~] = Circular_Object_Uncertainty(Tree_stem_geometry_init(1), Tree_stem_geometry_init(2), Tree_stem_geometry_init(3), Scanner_loc_cell, x_fd_points_cell, y_fd_points_cell, beam_divergence, sigma_range_device, max_incidence_angle, Point_Duplication);            
                laser_radius_cell = Laser_Beamwidth(beam_exit_diameter, beam_divergence, x_fd_points_cell, y_fd_points_cell, Scanner_loc_cell, Point_Duplication);

                % The initial expected distance is used to normalise the objective or as weights
                Distance_Plot = false;
                Distance_Print = false;
                
                [~, Expected_distance_init_cell] = Expected_Circle_Distance(distance_moment, Scanner_loc_cell, x_fd_points_cell, y_fd_points_cell, sigma_radial_cell, sigma_range_cell, range_bias, laser_radius_cell, Tree_stem_geometry_init(1), Tree_stem_geometry_init(2), Tree_stem_geometry_init(3), Point_Mirroring, Point_Duplication, Distance_Plot, Distance_Print);
                
                %--% Solution bounds %--%
                % Reasonable bounds for the optimisation process are based on the radius
                position_margin = bounds_margin * Tree_stem_geometry_init(3);      
                radius_margin   = 1 + bounds_margin;

                geometry_ub = Tree_stem_geometry_init .* [1, 1, radius_margin] + [position_margin, position_margin, 0];
                geometry_lb = Tree_stem_geometry_init .* [1, 1, 1/radius_margin] - [position_margin, position_margin, 0];

                % The initial estimate is normalised
                Tree_stem_geometry_init_n = (Tree_stem_geometry_init - geometry_lb) ./ (geometry_ub - geometry_lb);
                
                % As a result, the optimisation bounds are between 0 and 1
                Optim_LB = zeros(1, length(Tree_stem_geometry_init));
                Optim_UB = ones(1, length(Tree_stem_geometry_init));

                %--% Optimisation %--%
                % Optimisation of both likelihood and expected distance
                Tree_stem_geometry_steps_n  = [];       % The optimisation process is saved
                Objective_value_steps_n     = [];
                
                Tree_stem_optimisation = @(Tree_stem_geometry_n) Circle_Fitting_Optimiser(objective_balance, distance_moment, Tree_stem_geometry_n, geometry_lb, geometry_ub, x_fd_points_cell, y_fd_points_cell, Expected_distance_init_cell, sigma_radial_cell, sigma_range_cell, laser_radius_cell, Scanner_loc_cell, range_bias, Point_Duplication, Point_Mirroring, Weighting, central_theta_list, delta_inner_list, delta_outer_list);
                
                % The optimisation options
                Options                 = optimoptions('fmincon', 'OutputFcn', @Optim_Path, 'optimalitytolerance', 1e-6, 'steptolerance', 1e-6, 'Display', 'Off');
            
                % Optimisation
                [Tree_stem_geometry_n, Objective_value, ~] = fmincon(Tree_stem_optimisation, Tree_stem_geometry_init_n, [], [], [], [], Optim_LB, Optim_UB, [], Options);
                
                %--% Result evaluation %--%
                % Check if any of the bounds are reached
                if ~isempty(find(Tree_stem_geometry_n == 0, 1)) || ~isempty(find(Tree_stem_geometry_n == 1, 1))
                    disp('Warning: the bounds for optimisation have restricted the solution. Try increasing the bounds margin. The script will continue after a button press')
                    pause()
                end
                
                % The geometry is non-normalised
                Tree_stem_geometry_fit    = Tree_stem_geometry_n .* (geometry_ub - geometry_lb) + geometry_lb;

                % The optimisation steps are non-normalised
                Tree_stem_geometry_steps    = Tree_stem_geometry_steps_n .* (geometry_ub - geometry_lb) + geometry_lb;
                Objective_value_steps       = Objective_value_steps_n;      % Matlab doesn't like global variables as outputs

                % The percentage change in geometry w.r.t. the initial radius (such that it is irrespective of the magnitudes of the coordinates)
                geometry_change_list    = abs((Tree_stem_geometry_fit - Tree_stem_geometry_init)) ./ Tree_stem_geometry_init(3) * 100;
                                
                % Check for convergence
                max_change              = max(geometry_change_list);
                
                if max_change < convergence_threshold                   % The maximum change is within the threshold; convergence is achieved
                    convergence = true;

                elseif Objective_value > Objective_value_old            % The current fit is worse than the previous one, so it is not converging to the right solution
                    convergence             = true;

                    Tree_stem_geometry_fit  = Tree_stem_geometry_init;  % The previous geometry is used
                                                        
                    Objective_value         = Objective_value_old;
                    
                else                                                    % The initial guess will be updated based on the generated geometry
                    Tree_stem_geometry_init = Tree_stem_geometry_fit;
                    
                    Objective_value_old     = Objective_value;
                end
            end
            
            % The fitted tree centre and radius
            Tree_centre_fit     = Tree_stem_geometry_fit(1 : 2);   
            Tree_radius_fit     = Tree_stem_geometry_fit(3);

        elseif strcmp(Fitting_method(1:2), 'LS')
            %--% Fitting with regular or Pratt least-squares %--%
            [Tree_centre_x_fit, Tree_centre_y_fit, Tree_radius_fit] = Least_Squares_Circle_Fitting(x_fd_points_cell, y_fd_points_cell, Fitting_method, Weighting, Scanner_loc_cell, beam_divergence, sigma_range_device, range_bias, max_incidence_angle, Point_Duplication);
            Tree_centre_fit = [Tree_centre_x_fit, Tree_centre_y_fit];
            
            %--% The objective function value %--%
            % The uncertainty of the point cloud given this new geometry
            [sigma_radial_cell, sigma_range_cell, ~, ~] = Circular_Object_Uncertainty(Tree_centre_x_fit, Tree_centre_y_fit, Tree_radius_fit, Scanner_loc_cell, x_fd_points_cell, y_fd_points_cell, beam_divergence, sigma_range_device, max_incidence_angle, Point_Duplication);                    
            laser_radius_cell = Laser_Beamwidth(beam_exit_diameter, beam_divergence, x_fd_points_cell, y_fd_points_cell, Scanner_loc_cell, Point_Duplication);

            % The normalised likelihood
            Diagnostics_Likelihood = false;
            [~, ~, Avg_max_likelihood, ~, ~, ~] = Circle_Fit_Likelihood(x_fd_points_cell, y_fd_points_cell, sigma_radial_cell, sigma_range_cell, laser_radius_cell, Tree_centre_x_fit, Tree_centre_y_fit, Tree_radius_fit, Scanner_loc_cell, range_bias, Point_Duplication, Point_Mirroring, Diagnostics_Likelihood);              
            Norm_expected_distance = 1;         % Note that it is 1, as there is no initial estimate to normalise it with
            
            [Geometry_Function_Handle, Geometry_variable_bounds] = Shape_Basis_Function();
            Objective_value = Objective_Function(Avg_max_likelihood, Norm_expected_distance, objective_balance, Geometry_Function_Handle, [Tree_centre_fit, Tree_radius_fit], Geometry_variable_bounds, central_theta_list, delta_inner_list, delta_outer_list, Scanner_loc_cell);

            % There is no solution process, so the outputs are left empty
            Tree_stem_geometry_steps    = [];
            Objective_value_steps       = [];
        end
                
    %% Fitting assessment %%       
        % Printed statements
        if Print == true
            fprintf('The circle has been fitted with an objective value of %.4g \n', Objective_value);

            if Objective_value > 1e1
                disp('The high objective value indicates that constraints are violated by the solution')
            end
        end      
            
        % Plot showing the point cloud and fitted shape
        if Diagnostics == true
            % The fitted circle parameters
            Tree_centre_x_fit = Tree_centre_fit(1);
            Tree_centre_y_fit = Tree_centre_fit(2);
            
            disp('----------')
            fprintf('   Centre of the fitted circle: [%g, %g] m \n', Tree_centre_x_fit, Tree_centre_y_fit);
            fprintf('   Radius of the fitted circle: %g m \n', Tree_radius_fit);

            % Circle coordinates
            theta_list      = linspace(0, 2*pi, 1e3);

            x_fit_list      = Tree_radius_fit * sin(theta_list) + Tree_centre_x_fit;
            y_fit_list      = Tree_radius_fit * cos(theta_list) + Tree_centre_y_fit;

            % Aggregated point cloud coordinates
            x_points        = vertcat(x_points_cell{:});
            y_points        = vertcat(y_points_cell{:});
            x_f_points      = vertcat(x_f_points_cell{:});
            y_f_points      = vertcat(y_f_points_cell{:});  
            x_fd_points     = vertcat(x_fd_points_cell{:});
            y_fd_points     = vertcat(y_fd_points_cell{:});

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
            set(gcf, 'color', [1, 1, 1])

            % Circles and point clouds
            hold on
            grid on

            plot(x_fit_list, y_fit_list, 'color', 'k', 'DisplayName', 'Fitted circle', 'LineWidth', 2);

            scatter(x_points, y_points, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', 'Original point cloud', 'MarkerFaceAlpha', 1.0);

            if Filtering == true
                if Point_Duplication == false
                    scatter(x_f_points, y_f_points, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Filtered point cloud', 'MarkerFaceAlpha', 1.0);
                else
                    scatter(x_fd_points, y_fd_points, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Filtered and duplicated point cloud', 'MarkerFaceAlpha', 1.0); 
                end
            else
                if Point_Duplication == true
                    scatter(x_fd_points, y_fd_points, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Duplicated point cloud', 'MarkerFaceAlpha', 1.0);                         
                end
            end

            % Axis limits
            x_lim = Tree_centre_x_fit + 2 * Tree_radius_fit * [-1, 1];
            y_lim = Tree_centre_y_fit + 2 * Tree_radius_fit * [-1, 1];

            xlim(x_lim);
            ylim(y_lim);

            % The aspect ratio
            AR = (max(y_lim) - min(y_lim)) / (max(x_lim) - min(x_lim));
            pbaspect([1, AR, 1])

            % Axis looks
            xlabel('x [m]');
            ylabel('y [m]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            % Legend
            legend('show', 'Location', 'Eastoutside');
            hold off

            disp('The script will continue (and close the plot) if you press a key on the keyboard')
            pause();           

            close(1);            
        end
        
    %% Optimisation path function %%
        function stop = Optim_Path(geometry_estimate_n, optimisation_values, ~)
            stop = false;
            
            % Concatenate current point and objective function
            Tree_stem_geometry_steps_n  = [Tree_stem_geometry_steps_n; geometry_estimate_n];
            Objective_value_steps_n     = [Objective_value_steps_n; optimisation_values.fval];
        end
end
      