% The Monte Carlo loop is used to determine the uncertainty in the cylinder geometry
% It can be ran in parallel if Parallel_Pool is true

function [Cylinder_Geometry, Gaussian_Mixture_Models, Monte_Carlo_Loop] = Monte_Carlo_Loop_Cylinder(Parallel_Pool, Scanner_Parameters, Scanning_Parameters, Monte_Carlo_Inputs, Statistical_Values, Fitting_Parameters, Point_Cloud_Coord, Output_Decisions)

    %% Input structures %%
        % Parallel pool
        Parallel_Loop               = Parallel_Pool.Parallel_Loop;
        idle_timeout                = Parallel_Pool.idle_timeout;
    
        % Point cloud
        number_points_list          = Point_Cloud_Coord.number_points_list;
        point_cloud_cell            = Point_Cloud_Coord.point_cloud_cell;

        % Fitting parameters
        distance_moment             = Fitting_Parameters.distance_moment;

        % Scanning parameters
        Scanner_loc_cell            = Scanning_Parameters.Scanner_loc_cell;
        number_scanners             = Scanning_Parameters.number_scanners;
        Initial_Uncertainty         = Scanning_Parameters.Initial_Uncertainty;

        % Scanner parameters
        range_bias                  = Scanner_Parameters.range_bias;

        % Statistical values
        Confidence_interval         = Statistical_Values.Confidence_interval;
        Bhatt_coeff_threshold       = Statistical_Values.Bhatt_coeff_threshold;
        
        % Monte Carlo simulation
        max_MC_length               = Monte_Carlo_Inputs.max_MC_length;  
        number_GMM_fitment_iter     = Monte_Carlo_Inputs.number_GMM_fitment_iter;
        
        % Outputs
        Print                       = Output_Decisions.Print;
        Output_Decisions_MC         = Output_Decisions;
        
    %% Manual inputs %%
        % Least-squares fitting is used for the initial cylinder fit and its uncertainty
        LS_Plot                         = false;        % [true, false]
        LS_Diagnostics                  = false;        % [true, false]
        Distr_Diagnostics               = false;        % [true, false]

        % Plots and progress messages are generally not desired within the MC loop
        Sampling_Diagnostics            = false;        % [true, false]
        Output_Decisions_MC.Print       = false;        % [true, false]
        Output_Decisions_MC.Plot        = false;        % [true, false]
        
        % Diagnostic information regarding the Monte Carlo loop's results
        Diagnostics_Mahal               = false;        % [true, false] 
        GMM_Diagnostics                 = false;        % [true, false]
        Monte_Carlo_Diagnostics         = false;        % [true, false]

        % Iterative search is used to determine the confidence intervals
        max_F_error                     = 1e-3;         % [-] Convergence threshold for the cumulative density function
        max_iterations                  = 1e2;          % [-] Due to convexity of the function convergence should be rapid, but this sets a maximum
        Search_Print                    = false;        % [true, false] Prints possible warnings and the outcome of the search

        % The top and bottom heights are assumed independent and a warning is printed if their correlation is greater than this threshold
        height_correlation_threshold    = 0.25;         % [-]
        
    %% Centroid of the point cloud %%
        % To improve robustness and sensibility of the results, geometry is relative to the original point cloud centroid
        point_cloud_matrix      = vertcat(point_cloud_cell{:});
        point_cloud_centroid    = mean(point_cloud_matrix, 1);
        num_dim                 = length(point_cloud_centroid);
        
    %% Initial uncertainty estimate %%
        % In order to sample the point cloud, its uncertainty must be estimated which is a product of the geometry
        % The initial geometry is determined through least-squares using a matrix of the point cloud
        LS_Geometry = Least_Squares_Cylinder_Fitting(Point_Cloud_Coord, Scanning_Parameters, Fitting_Parameters, Confidence_interval, LS_Diagnostics, LS_Diagnostics);

        [cylinder_centre_init, cylinder_radius_init, cylinder_direction_init, cylinder_length_init] = deal(LS_Geometry.Cylinder_centre, LS_Geometry.Cylinder_radius, LS_Geometry.Cylinder_direction, LS_Geometry.Cylinder_length);

        % To robustly reduce the dimension of the cylinder centre through projection, the largest component is noted of the initial vector
        max_component_boolean   = false(1, num_dim);
        max_component           = max(abs(cylinder_direction_init));
        max_component_ind       = find(abs(cylinder_direction_init) == max_component, 1);       % Ensures that only 1 component is chosen if multiple are equal in magnitude   
        max_component_boolean(max_component_ind) = true;
        
        % The uncertainty of each point, given the geometry
        [sigma_radial_cell, sigma_prop_cell, ~, ~] = Cylindrical_Object_Uncertainty(cylinder_centre_init, cylinder_direction_init, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);

        alpha                       = 1 - Confidence_interval/100;
        Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

        % Plot to check the geometry and uncertainty derived from it
        if LS_Plot == true
            % The point cloud distributions
            Diagnostics_Ellipsoids      = false;            
            Point_Cloud_Distributions_c = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Diagnostics_Ellipsoids);

            % Colour map for the points of each scanner
            scanner_cmap    = cbrewer('qual', 'Set2', max(3, number_scanners));
            
            % The number of coordinates in the cylinder
            number_coord = 1e3;     
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
                
            % Initial least-squares cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius_init, cylinder_length_init, cylinder_centre_init, cylinder_direction_init, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial least-squares cylinder');

            % The point cloud and its uncertainty
            distribution_mu_cell    = Point_Cloud_Distributions_c.distribution_mu_cell;
            distribution_radii_cell = Point_Cloud_Distributions_c.distribution_radii_cell;
            distribution_axes_cell  = Point_Cloud_Distributions_c.distribution_axes_cell;
            num_distributions       = Point_Cloud_Distributions_c.num_distributions;
            
            number_beams_list       = Point_Cloud_Coord.number_beams_list;
            number_cumulative_beams = [0, cumsum(number_beams_list)];
            
            for d = 1: num_distributions
                % This distribution's scanner and its colour
                scanner_ind     = find(number_cumulative_beams - d >= 0, 1) - 1;   
                scanner_colour  = scanner_cmap(scanner_ind, :);
                
                % Distribution properties
                distribution_mu     = distribution_mu_cell{d};      
                distribution_radii  = distribution_radii_cell{d};
                distribution_axes   = distribution_axes_cell{d};
                
                % Surface
                [ellipsoid_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_coord);

                x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', scanner_colour, 'FaceAlpha', 0.10, 'DisplayName', sprintf('Point cloud uncertainty, \\alpha = %.3g', alpha));

                % Centre (mu)
                sc_mu = scatter3(distribution_mu(1), distribution_mu(2), distribution_mu(3), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', scanner_colour, 'DisplayName', '\mu');

                if d > 1
                    surf_el.HandleVisibility    = 'Off';
                    sc_mu.HandleVisibility      = 'Off';
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

            disp('The script will continue and the plot will close when a key is pressed');
            pause();

            close(1);             
        end
        
    %% Monte Carlo loop %%
        % Initiate the parallel pool
        if Parallel_Loop == true
            % The number of maximum iterations is a multiple of the number of cores
            number_cores    = feature('numcores');                  % Checks the number of available cores
            number_loops    = ceil(max_MC_length / number_cores);
            max_MC_length   = number_cores * number_loops;    
            
            % The parallel pool is started
            Parallel_Pool_Starter(idle_timeout, number_cores)
            number_cores_parallel = number_cores;
        else
            % Otherwise the number of cores used in parallel are set to 0 to make the parfor run as a for
            number_cores_parallel   = 0;
            number_cores            = 1;
            number_loops            = max_MC_length;
        end
        
        % As the data is resampled, the uncertainty is approximately twice what it would normally be if there was initial uncertainty and the scanner specifications are changed accordingly
        Scanner_Parameters_MC = Scanner_Parameters;

        if Initial_Uncertainty == true
            Scanner_Parameters_MC.beam_divergence       = 2*Scanner_Parameters.beam_divergence;
            Scanner_Parameters_MC.beam_exit_diameter    = 2*Scanner_Parameters.beam_exit_diameter;
            Scanner_Parameters_MC.sigma_range_device    = 2*Scanner_Parameters.sigma_range_device;
        end

        % The parameters are saved every iteration for both the initial least-squares geometry and the final optimal geometry
        MC_Opt_Geometry_Data        = struct('Cylinder_centre', [], 'Cylinder_direction', [], 'Cylinder_radius', [], 'Cylinder_length', [], 'Cylinder_height_top', [], 'Cylinder_height_bot', []);
        MC_LS_Geometry_Data         = struct('Cylinder_centre', [], 'Cylinder_direction', [], 'Cylinder_radius', [], 'Cylinder_length', [], 'Cylinder_height_top', [], 'Cylinder_height_bot', []);
        
        % The expected Mahalanobis distance, incidence angle, propagation error and cylinder distance for each point are saved each iteration
        MC_expected_Mahal_dist_cell = cell(1, max_MC_length);
        MC_incidence_angle_cell     = cell(1, max_MC_length);
        MC_propagation_error_cell   = cell(1, max_MC_length);
        MC_cylinder_distance_cell   = cell(1, max_MC_length);

        % The optimiser diagnostics are saved as well for the infinite cylinder fitting
        MC_Optimiser_Diagnostics    = struct('Optimiser_Steps', [], 'Optimum_First_Derivatives', [], 'Norm_Opt_First_Derivatives', [], 'Optimum_Second_Derivatives', [], 'Norm_Opt_Second_Derivatives', [], 'Optimum_Hessian', [], 'Norm_Optimum_Hessian', [], 'parameter_labels', [], 'parameter_units', [], 'num_parameters', [], 'number_UU_iterations', []);

        % To later determine properties of the models, they are saved in a structure
        MC_Gaussian_Mixture_Models  = struct('max_component_boolean', max_component_boolean, 'GMM', [], 'Data', struct('Proj_cylinder_centre', [], 'Cylinder_direction_xy', [], 'Cylinder_radius', [], 'Cylinder_height_top', [], 'Cylinder_height_bot', []));

        % Progress counter (to maximum)
        DQ      = parallel.pool.DataQueue;
        tick    = 0;
        N       = max_MC_length;
        afterEach(DQ, @ProgressUpdate);
    
        Convergence     = false;
        loop_counter    = 0;

        while Convergence == false && loop_counter < number_loops
            % The loop counter is updated until convergence or the maximum is reached
            loop_counter    = loop_counter + 1;
            
            start_ind       = (loop_counter - 1) * number_cores + 1;
            end_ind         = loop_counter * number_cores;
            
            % Parallel loop            
            GMM_Variable_Data = MC_Gaussian_Mixture_Models.Data;              % Temporary structure which can be used in the parallel loop

            parfor (c = start_ind : end_ind, number_cores_parallel)
                %--% Point cloud sampling %--%
                % Set the RNG seed for consistency (in a parallel loop the initially specified RNG seed is ignored)
                rng(c + 100);

                % Sample the uncertain point cloud
                Point_Cloud_Coord_MC    = Gaussian_Point_Sampling(Point_Cloud_Distributions, Sampling_Diagnostics);
                
                point_cloud_matrix_MC   = vertcat(Point_Cloud_Coord_MC.point_cloud_cell{:});
                point_cloud_centroid_MC = mean(point_cloud_matrix_MC, 1);

                %--% Cylinder fitting %--%
                [Optimal_Geometry, LS_Geometry, Optimiser_Diagnostics, ~] = Cylinder_Geometry_Estimation(Point_Cloud_Coord_MC, Scanning_Parameters, Scanner_Parameters_MC, Statistical_Values, Fitting_Parameters, Output_Decisions_MC);
                
                [cylinder_centre, cylinder_direction, cylinder_radius, cylinder_top_loc, cylinder_bot_loc, cylinder_length] = deal(Optimal_Geometry.Cylinder_centre, Optimal_Geometry.Cylinder_direction, Optimal_Geometry.Cylinder_radius, Optimal_Geometry.Top_loc, Optimal_Geometry.Bottom_loc, Optimal_Geometry.Cylinder_length);
                
                % The top and bottom locations are converted to heights
                cylinder_height_top     = (cylinder_top_loc - point_cloud_centroid) / cylinder_direction;          % Note that the height is relative to the original point cloud centroid
                cylinder_height_bot     = (cylinder_bot_loc - point_cloud_centroid) / cylinder_direction;

                cylinder_direction_LS   = LS_Geometry.Cylinder_direction;
                cylinder_height_list_LS = ([LS_Geometry.Top_loc; LS_Geometry.Bottom_loc] - point_cloud_centroid) / cylinder_direction_LS;      
                cylinder_height_bot_LS  = min(cylinder_height_list_LS);                     % To ensure consistency
                cylinder_height_top_LS  = max(cylinder_height_list_LS);
                
                % The heights are appended to the structures, and the old locations are removed as well as other unneeded fields
                Optimal_Geometry.Cylinder_height_top    = cylinder_height_top;
                Optimal_Geometry.Cylinder_height_bot    = cylinder_height_bot;
                Optimal_Geometry                        = rmfield(Optimal_Geometry, {'Top_loc', 'Bottom_loc'});
                
                LS_Geometry.Cylinder_height_top         = cylinder_height_top_LS;
                LS_Geometry.Cylinder_height_bot         = cylinder_height_bot_LS;
                LS_Geometry                             = rmfield(LS_Geometry, {'Top_loc', 'Bottom_loc'});
                
                % The cylinder centre of the optimiser is defined w.r.t. its sampled point cloud centroid, here it is defined w.r.t. the original centroid 
                cylinder_axis_start     = cylinder_centre - point_cloud_centroid_MC + point_cloud_centroid;
                [cylinder_centre, ~, ~] = Point_to_Vector_Projection(point_cloud_centroid, cylinder_direction, cylinder_axis_start);

                % The vectors are made consistent by ensuring the 3rd component is positive
                cylinder_direction                      = sign(cylinder_direction(num_dim)) * cylinder_direction;
                Optimal_Geometry.Cylinder_direction     = cylinder_direction;
                
                cylinder_direction_LS                   = sign(cylinder_direction_LS(num_dim)) * cylinder_direction_LS;
                LS_Geometry.Cylinder_direction          = cylinder_direction_LS;
                
                % The results are saved for later analysis
                MC_Opt_Geometry_Data(c)     = Optimal_Geometry;
                MC_LS_Geometry_Data(c)      = LS_Geometry;
                MC_Optimiser_Diagnostics(c) = Optimiser_Diagnostics;

                %--% Expected Mahalanobis distance for each point %--%
                [MC_sigma_rad_cell, MC_sigma_prop_cell, ~, ~]   = Cylindrical_Object_Uncertainty(cylinder_centre, cylinder_direction, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);
                MC_Point_Cloud_Distributions                    = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, MC_sigma_rad_cell, MC_sigma_prop_cell, Scanner_loc_cell, range_bias, Diagnostics_Mahal);

                [~, expected_Mahal_dist_list]   = Expected_Mahal_Distance_Cylinder_Line_Approx(cylinder_centre, cylinder_radius, cylinder_direction, distance_moment, MC_Point_Cloud_Distributions, Scanner_loc_cell, Diagnostics_Mahal);
                MC_expected_Mahal_dist_cell{c}  = expected_Mahal_dist_list;
            
                %--% Propagation error for each point %--%
                % The scanner locations are repeated for their point cloud
                number_points_cell  = num2cell(number_points_list);
                Scanner_Rep_fun     = @(number_points, scanner_loc) repmat(scanner_loc, [number_points, 1]);
                scanner_loc_cell    = cellfun(Scanner_Rep_fun, number_points_cell, Scanner_loc_cell, 'UniformOutput', false);
                scanner_loc_matrix  = vertcat(scanner_loc_cell{:});

                % Vectors from the scanners to each point
                scan_vector_matrix  = point_cloud_matrix - scanner_loc_matrix;
                point_range_list    = sqrt(sum(scan_vector_matrix.^2, 2));
                scan_vector_matrix  = scan_vector_matrix ./ point_range_list;

                % Intersections for each scanning vector with the cylinder
                [intersection_point_matrix, cylinder_intersects, ~, ~, incidence_angle_list]    = Vector_Cylinder_Intersection(scanner_loc_matrix, scan_vector_matrix, cylinder_radius, cylinder_centre, cylinder_length, cylinder_direction);
                MC_incidence_angle_cell{c}                                                      = incidence_angle_list;

                % The propagation errors. Non-intersects are given a NaN value
                intersect_vector_matrix = intersection_point_matrix - scanner_loc_matrix;
                intersect_range_list    = sqrt(sum(intersect_vector_matrix.^2, 2));

                propagation_error_list                          = point_range_list - intersect_range_list;
                propagation_error_list(~cylinder_intersects)    = NaN;

                MC_propagation_error_cell{c}                    = propagation_error_list;

                %--% Signed cylinder distance for each point %--%
                % The points are projected onto the cylinder axis
                [~, ~, omega_list]              = Point_to_Vector_Projection(point_cloud_matrix_MC, cylinder_direction, cylinder_centre);
                cylinder_distance_list          = omega_list - cylinder_radius;         % Note that a point lying outside the cylinder is positive, inside negative

                MC_cylinder_distance_cell{c}    = cylinder_distance_list;

                %--% Gaussian mixture model data %--%     
                % The first two components are used, as the third one is always positive  
                cylinder_direction_xy   = cylinder_direction(1 : num_dim - 1);
                
                % The projected 2D circle centre
                cylinder_direction_k    = cylinder_direction(max_component_boolean);                                        % Largest component
                cylinder_direction_ij   = cylinder_direction(~max_component_boolean);                                       % Remaining two components
                
                cylinder_centre_c       = cylinder_centre - point_cloud_centroid;                                           % To improve robustness, the centroid of the point cloud is removed
                lambda                  = cylinder_centre_c(max_component_boolean) / cylinder_direction_k;                  % Distance along vector to the plane of the two smallest vector components
                proj_cylinder_centre    = cylinder_centre_c(~max_component_boolean) - lambda * cylinder_direction_ij;       % 2D cylinder centre projected on said plane

                GMM_Variable_Data(c) = struct('Proj_cylinder_centre', proj_cylinder_centre, 'Cylinder_direction_xy', cylinder_direction_xy, 'Cylinder_radius', cylinder_radius, 'Cylinder_height_top', cylinder_height_top, 'Cylinder_height_bot', cylinder_height_bot);
                
                % Progress update
                send(DQ, c);
            end

            % The original structure is changed according to the temporary one
            MC_Gaussian_Mixture_Models.Data = GMM_Variable_Data;

            %--% Gaussian mixture model fitment %--%
            % The geometry data thus far
            [current_geometry_matrix, ~] = Structure_Data_Concatenation(MC_Gaussian_Mixture_Models.Data);

            % The optimal models are fitted based on the AIC
            [GM_Model, Shared_Covariance, number_GM_components, AICc_min] = Gaussian_Mixture_Model_Fitting(current_geometry_matrix, number_GMM_fitment_iter);
            MC_Gaussian_Mixture_Models.GMM{loop_counter} = GM_Model;

            %--% Bhattacharyya coefficient for convergence %--%
            if loop_counter > 1 && ~isempty(GM_Model)        % If a GMM could not be fitted, empty values are given
                % The previous loop's Gaussian mixture model
                GMM_prev = MC_Gaussian_Mixture_Models.GMM{loop_counter - 1};

                % The Bhattacharyya coefficient is computed between the current and previous models
                if ~isempty(GMM_prev)
                    [Bhatt_coeff, Bhatt_distance] = Bhattacharyya_GMM(GM_Model, GMM_prev);

                % If a previous model does not exist, NaN values are given
                else
                    [Bhatt_coeff, Bhatt_distance] = deal(NaN);
                end
            else
                [Bhatt_coeff, Bhatt_distance] = deal(NaN);
            end

            % All the GMM metrics are saved
            MC_Gaussian_Mixture_Models.GMM_Metrics(loop_counter) = struct('Shared_Covariance', Shared_Covariance, 'number_GM_components', number_GM_components, 'AICc_min', AICc_min, 'Bhatt_coeff', Bhatt_coeff, 'Bhatt_distance', Bhatt_distance);
                                          
            %--% Convergence check %--%            
            % It has converged if the coefficient is above the threshold            
            number_MC_iterations = end_ind;

            if Bhatt_coeff > Bhatt_coeff_threshold
                Convergence = true;
                
                % Convergence information is printed if desired
                if Print == true
                    fprintf('\n');
                    disp('---------------------------');
                    fprintf('%g Monte Carlo iterations were performed \n', number_MC_iterations);
                    fprintf('The GMM has a Bhattacharyya coefficient of %.3g \n', Bhatt_coeff);
                    disp('---------------------------');
                    fprintf('\n');
                end
            end
        end

        % The expected Mahalanobis distance, incidence angle, propagation error and cylinder distances are aggregated into matrices
        MC_expected_Mahal_dist_matrix   = horzcat(MC_expected_Mahal_dist_cell{:});
        MC_incidence_angle_matrix       = horzcat(MC_incidence_angle_cell{:});
        MC_propagation_error_matrix     = horzcat(MC_propagation_error_cell{:});
        MC_cylinder_distance_matrix     = horzcat(MC_cylinder_distance_cell{:});

        % A warning or error is displayed if the model has not converged within the number of steps
        if Convergence == false
            if isnan(Bhatt_coeff)
                error('The Gaussian mixture model could not be created within %g iterations. \n', max_MC_length);
            else
                warning('The algorithm could not converge within %g iterations. The Bhattacharyya coefficient is still %.3g \n', max_MC_length, Bhatt_coeff);
            end            
        end
        
        % The Monte Carlo data is collected in a single structure
        Monte_Carlo_Loop = struct('Convergence', Convergence, 'number_MC_iterations', number_MC_iterations, 'MC_Geometry_Data', struct('Optimal', MC_Opt_Geometry_Data, 'LS', MC_LS_Geometry_Data), 'Optimiser_Diagnostics', MC_Optimiser_Diagnostics, 'Gaussian_Mixture_Models', MC_Gaussian_Mixture_Models, 'expected_Mahal_dist_matrix', MC_expected_Mahal_dist_matrix, 'incidence_angle_matrix', MC_incidence_angle_matrix, 'propagation_error_matrix', MC_propagation_error_matrix, 'cylinder_distance_matrix', MC_cylinder_distance_matrix);

    %% Final Gaussian mixture model %%
        % The correlation and covariance between all variables are computed for the method of moments
        [~, GMM_Variable_Correlation] = Correlation_GMM(GM_Model, GMM_Variable_Data);
    
        % 1D Gaussian distributions for the expected values and uncertainties
        try
            [GMM_Variable_Distributions, ~, ~, ~] = Independent_Gaussians_GMM(GM_Model, GMM_Variable_Data, Confidence_interval, max_F_error, max_iterations, Search_Print, GMM_Diagnostics);
        catch            
            % If a multivariate-normal matrix is not positive-definite, analytical projection fails and instead the values are computed numerically
            warning('The GMM contains an MVN matrix that is not positive-definite. Numerical projection is performed instead.');
            GMM_Variable_Distributions = Independent_Gaussians_GMM_Numerical(GMM_Variable_Data, Confidence_interval);
        end
                
        % Complete structure
        Gaussian_Mixture_Models = struct('max_component_boolean', max_component_boolean, 'GMM', GM_Model, 'GMM_Metrics', MC_Gaussian_Mixture_Models.GMM_Metrics(loop_counter), 'Distributions', GMM_Variable_Distributions, 'Correlation', GMM_Variable_Correlation, 'Data', GMM_Variable_Data);
    
    %% Monte Carlo loop diagnostics %%  
        if Monte_Carlo_Diagnostics == true              
            %--% GMM metrics %--%
            metric_names    = fieldnames(MC_Gaussian_Mixture_Models(1).GMM_Metrics);
            number_metrics  = length(metric_names);
            
            metric_colours  = cbrewer('qual', 'Set1', max(number_metrics, 3));    
            metric_colours  = max(metric_colours, 0);
            metric_colours  = min(metric_colours, 1);
                            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    

            hold on

            for m = 1 : number_metrics
                % This metric's data
                metric_name     = metric_names{m};
                metric_colour   = metric_colours(m, :);

                metric_cell     = {MC_Gaussian_Mixture_Models.GMM_Metrics.(metric_name)};
                metric_data     = horzcat(metric_cell{:});

                % Its subplot
                subplot(2, ceil(number_metrics / 2), m)
                grid on
                hold on

                scatter(1 : loop_counter, metric_data, 'filled', 'markerfacecolor', metric_colour);

                xlim([1, loop_counter]);
                xlabel('Loop [-]');

                metric_str = strrep(metric_name, '_', ' ');
                ylabel(sprintf('%s [-]', metric_str));

                % Font size
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            end
            
            %--% Cumulative distribution functions %--%
            % The data used to form the Gaussian mixture model
            [data_matrix, ~] = Structure_Data_Concatenation(MC_Gaussian_Mixture_Models.Data);
                
            % The CDF values for each loop
            GMM_CDF_cell = cell(1, loop_counter);
                        
            for n = loop_counter : -1 : 1
                % The Gaussian mixture model
                GM_Model = MC_Gaussian_Mixture_Models.GMM{n};

                if isempty(GM_Model)
                    continue
                else
                    % Their cumulative density values
                    GMM_CDF_list = cdf(GM_Model, data_matrix);

                    % Sorted based on the last CDF values
                    if n == loop_counter
                        [~, order]  = sort(GMM_CDF_list);
                    end

                    GMM_CDF_cell{n} = GMM_CDF_list(order);
                end
            end
            
            % Colours used for the GMM CDFs
            CDF_colour_matrix = cbrewer('seq', 'OrRd', max(loop_counter, 3));
            CDF_colour_matrix = max(0, CDF_colour_matrix);
            CDF_colour_matrix = min(1, CDF_colour_matrix);

            figure(2)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    

            hold on

            sgtitle('Later means darker');
                
            % The first subplot contains only the last two GMMs' CDF values (i.e. the ones that are supposed to be similar), the second contains all
            loop_indices_cell   = {loop_counter : -1 : loop_counter - 1, 1 : loop_counter};
            line_width_list     = [2, 1];
            number_subplots     = length(loop_indices_cell);
                
            for sp = 1 : number_subplots
                % This subplot's data
                loop_indices    = loop_indices_cell{sp};
                line_width      = line_width_list(sp);

                subplot(1, number_subplots, sp)
                hold on
                grid on

                for n = loop_indices
                    % This loop's data
                    GMM_CDF_list    = GMM_CDF_cell{n};
                    CDF_colour      = CDF_colour_matrix(n, :);

                    if ~isempty(GMM_CDF_list)
                        plot(GMM_CDF_list, 'color', CDF_colour, 'LineWidth', line_width);                    
                    end
                end

                % Axis labels
                xlabel('Sorted sample [-]');
                ylabel('F(X) [-]');

                % Font size
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            end
            
            %--% First and second derivatives at optimum %--%
            Optimiser_Diagnostics_cell  = num2cell(Monte_Carlo_Loop.Optimiser_Diagnostics);
                
            num_optim_parameters    = Optimiser_Diagnostics_cell{1}.num_parameters;
            optim_parameter_labels  = Optimiser_Diagnostics_cell{1}.parameter_labels;
            optim_parameter_units   = Optimiser_Diagnostics_cell{1}.parameter_units;

            parameter_colours    = cbrewer('qual', 'Set1', max(num_optim_parameters, 3));
            parameter_colours    = max(parameter_colours, 0);
            parameter_colours    = min(parameter_colours, 1);

            First_Derivatives_fun           = @(Optimiser_Diagnostics) Optimiser_Diagnostics.Optimum_First_Derivatives.first_derivative_list;
            optim_first_derivatives_cell    = cellfun(First_Derivatives_fun, Optimiser_Diagnostics_cell, 'UniformOutput', false);
            optim_first_derivatives_matrix  = vertcat(optim_first_derivatives_cell{:});

            Second_Derivatives_fun          = @(Optimiser_Diagnostics) Optimiser_Diagnostics.Optimum_Second_Derivatives.second_derivative_list;
            optim_second_derivatives_cell   = cellfun(Second_Derivatives_fun, Optimiser_Diagnostics_cell, 'UniformOutput', false);
            optim_second_derivatives_matrix = vertcat(optim_second_derivatives_cell{:});

            figure(3)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    

            for p = 1 : num_optim_parameters
                % This parameter's data
                parameter_colour        = parameter_colours(p, :);
                parameter_label         = optim_parameter_labels{p};
                parameter_unit          = optim_parameter_units{p};
                
                parameter_first_derivative_list = optim_first_derivatives_matrix(:, p);
                parameter_sec_derivative_list   = optim_second_derivatives_matrix(:, p);

                % First derivative histogram
                subplot(2, num_optim_parameters, p)
                hold on
                grid on

                [number_bins, ~] = Histogram_Bins(parameter_first_derivative_list);

                hist_grad   = histogram(parameter_first_derivative_list, number_bins, 'FaceAlpha', 1.0, 'FaceColor', parameter_colour, 'normalization', 'pdf');
                f_max       = max(hist_grad.Values);       % The maximum density value

                % Axes
                xlabel(sprintf('dQ/d%s [1/%s]', parameter_label, parameter_unit));              
                ylabel(sprintf('Probability density [1/%s]', parameter_unit));
                ylim([0, f_max]);
 
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                % Second derivative histogram
                subplot(2, num_optim_parameters, p + num_optim_parameters)
                hold on
                grid on

                [number_bins, ~] = Histogram_Bins(parameter_sec_derivative_list);

                hist_grad   = histogram(parameter_sec_derivative_list, number_bins, 'FaceAlpha', 1.0, 'FaceColor', parameter_colour, 'normalization', 'pdf');
                f_max       = max(hist_grad.Values);       % The maximum density value

                % Axes
                xlabel(sprintf('d^2Q/d%s^2 [1/%s^2]', parameter_label, parameter_unit));              
                ylabel(sprintf('Probability density [1/%s]', parameter_unit));
                ylim([0, f_max]);
 
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
            end

            %--% Pause %--%
            disp('The script will continue and the figures will close when a key is pressed.');
            pause()
            close(1 : 3);
        end
        
    %% Geometry parameters and their independent uncertainty %%
        % The expected values and uncertainty are saved in a structure, so that they can be interpreted more easily
        Cylinder_Geometry = struct();

        %--% Cylinder radius %--%
        % The radius properties can be copied over directly
        Cylinder_Geometry.Cylinder_radius = GMM_Variable_Distributions.Cylinder_radius;
        
        %--% Cylinder length %--%
        % As the top and bottom height are assumed independent, the expected value and confidence interval follow directly from their respective values
        height_correlation = GMM_Variable_Correlation.Cylinder_height_top.Cylinder_height_bot.correlation_coeff_matrix;
        
        if abs(height_correlation) > height_correlation_threshold
            warning('The top and bottom height are not uncorrelated, and their correlation coefficient is %.3g', height_correlation);
        end
        
        [height_top_mu, height_bot_mu]          = deal(GMM_Variable_Distributions.Cylinder_height_top.mu, GMM_Variable_Distributions.Cylinder_height_bot.mu);
        [height_top_sigma, height_bot_sigma]    = deal(GMM_Variable_Distributions.Cylinder_height_top.sigma, GMM_Variable_Distributions.Cylinder_height_bot.sigma);
    
        length_mu       = height_top_mu - height_bot_mu;
        length_sigma    = sqrt(height_top_sigma^2 + height_bot_sigma^2);            % Quadrature follows from the sum of independent variables
        
        m_STD                       = sqrt(chi2inv(Confidence_interval / 100, 1));
        length_confidence_interval  = length_mu + m_STD * [1; -1] * length_sigma;
        length_confidence_interval  = max(0, length_confidence_interval);           % As the length cannot be negative

        Cylinder_Geometry.Cylinder_length = struct('mu', length_mu, 'sigma', length_sigma, 'confidence_interval', length_confidence_interval);
    
        %--% Cylinder axis %--%
        % Due to the vector having a unit norm and the third component being defined as always positive, it can be computed from the other two vector components
        vector_xy_mu            = GMM_Variable_Distributions.Cylinder_direction_xy.mu;
        vector_z_mu             = sqrt(1 - sum(vector_xy_mu.^2));
        
        cylinder_direction_mu   = [vector_xy_mu, vector_z_mu];
        
        % They cannot be assumed independent, so the correlation coefficient is required for the uncertainty
        vector_xy_sigma = GMM_Variable_Distributions.Cylinder_direction_xy.sigma;
        rho_vector_xy   = GMM_Variable_Correlation.Cylinder_direction_xy.Cylinder_direction_xy.correlation_coeff_matrix(1, 2);               % The ij'th coefficient is off-diagonal
        
        vector_z_sigma  = 1/vector_z_mu * sqrt(sum((vector_xy_sigma.*vector_xy_mu).^2) + 2*rho_vector_xy*prod(vector_xy_mu)*prod(vector_xy_sigma));
        
        cylinder_direction_sigma = [vector_xy_sigma, vector_z_sigma];
        
        % The confidence interval values are ensured to be within possible bounds 
        vector_xy_confidence_interval   = GMM_Variable_Distributions.Cylinder_direction_xy.confidence_interval;
        vector_z_confidence_interval    = vector_z_mu + m_STD * [1; -1] * vector_z_sigma;
        
        cylinder_direction_confidence_interval = [vector_xy_confidence_interval, vector_z_confidence_interval];
        cylinder_direction_confidence_interval = min(cylinder_direction_confidence_interval, 1);
        cylinder_direction_confidence_interval = max(cylinder_direction_confidence_interval, 0);
        
        Cylinder_Geometry.Cylinder_direction = struct('mu', cylinder_direction_mu, 'sigma', cylinder_direction_sigma, 'confidence_interval', cylinder_direction_confidence_interval);

        %--% Cylinder centre %--%
        % The cylinder centre requires the projected centre, vector and height parameters and is thus not as simple to compute properties for
        % First the projected cylinder centre is converted to the original coordinate frame
        proj_cylinder_centre_mu                             = zeros(1, num_dim);
        proj_cylinder_centre_mu(~max_component_boolean)     = GMM_Variable_Distributions.Proj_cylinder_centre.mu;
              
        % Symbolic variables for the cylinder parameters used for the Gaussian mixture models
        % Note that the 3-dimensional projected cylinder centre is used as the zero component is removed later
        [proj_centre_x, proj_centre_y, proj_centre_z, vector_x, vector_y, x_bar, y_bar, z_bar, h_t, h_b] = deal([]);                    %#ok<ASGLU> Have to be initiated beforehand as this is a nested function
        syms proj_centre_x proj_centre_y proj_centre_z vector_x vector_y x_bar y_bar z_bar h_t h_b
    
        % These are substituted for the true values using the following function
        origin              = zeros(1, num_dim);        % As the point cloud has been centered
        eq_substitution_fun = @(eq) double(subs(eq, {proj_centre_x, proj_centre_y, proj_centre_z, vector_x, vector_y, x_bar, y_bar, z_bar, h_t, h_b}, num2cell([proj_cylinder_centre_mu, cylinder_direction_mu(1 : num_dim - 1), origin, height_top_mu, height_bot_mu])));
    
        % Expected value for the cylinder centre along the vector
        eq_cyl_centre_x = proj_centre_x + (h_t + h_b)/2 * vector_x + ((x_bar - proj_centre_x)*vector_x + (y_bar - proj_centre_y)*vector_y + (z_bar - proj_centre_z)*sqrt(1 - vector_x^2 - vector_y^2)) * vector_x;
        eq_cyl_centre_y = proj_centre_y + (h_t + h_b)/2 * vector_y + ((x_bar - proj_centre_x)*vector_x + (y_bar - proj_centre_y)*vector_y + (z_bar - proj_centre_z)*sqrt(1 - vector_x^2 - vector_y^2)) * vector_y;
        eq_cyl_centre_z = proj_centre_z + (h_t + h_b)/2 * sqrt(1 - vector_x^2 - vector_y^2) + ((x_bar - proj_centre_x)*vector_x + (y_bar - proj_centre_y)*vector_y + (z_bar - proj_centre_z)*sqrt(1 - vector_x^2 - vector_y^2)) * sqrt(1 - vector_x^2 - vector_y^2);
    
        cylinder_centre_mu_c    = cellfun(eq_substitution_fun, {eq_cyl_centre_x, eq_cyl_centre_y, eq_cyl_centre_z});
        cylinder_centre_mu      = cylinder_centre_mu_c + point_cloud_centroid;
 
        % To propagate uncertainty, the partial derivatives of the cylinder centre w.r.t. the other variables are computed
        Jacobian_eq_matrix  = jacobian([eq_cyl_centre_x; eq_cyl_centre_y; eq_cyl_centre_z], [proj_centre_x, proj_centre_y, proj_centre_z, vector_x, vector_y, h_t, h_b]);
        Jacobian_matrix     = arrayfun(eq_substitution_fun, Jacobian_eq_matrix);
    
        % Entries related to the zero-component of the projected cylinder centre are removed as its propagated uncertainty is zero
        Jacobian_matrix(:, [max_component_boolean, false(1, 3)]) = [];                            
        
        % The Jacobian is converted to structure format so that the parameter names can be used
        Jacobian_struct = struct('Proj_cylinder_centre', Jacobian_matrix(:, 1:2), 'Cylinder_direction_xy', Jacobian_matrix(:, 3:4), 'Cylinder_height_top', Jacobian_matrix(:, 5), 'Cylinder_height_bot', Jacobian_matrix(:, 6));
        
        % The uncertainty of each variable is propagated to the dimensions of the cylinder centre
        propagation_variables   = fieldnames(Jacobian_struct);
        number_prop_variables   = length(propagation_variables);
        
        uncert_prop_matrix      = zeros(number_prop_variables, number_prop_variables, num_dim);
        
        for i = 1 : number_prop_variables
            % The i-th variable's properties
            variable_i          = propagation_variables{i};
            sigma_list_i        = GMM_Variable_Distributions.(variable_i).sigma;        % Note that the variables may be multi-dimensional
            Jacobian_matrix_i   = Jacobian_struct.(variable_i);

            for j = i : number_prop_variables
                % The j-th variable's properties
                variable_j          = propagation_variables{j};
                sigma_list_j        = GMM_Variable_Distributions.(variable_j).sigma;
                Jacobian_matrix_j   = Jacobian_struct.(variable_j);

                % Their correlation coefficients matrix
                correlation_coeff_matrix_ij = GMM_Variable_Correlation.(variable_i).(variable_j).correlation_coeff_matrix;
                
                % Cycle through the dimensions of the cylinder centre
                for d = 1 : num_dim
                    part_der_list_i = Jacobian_matrix_i(d, :);
                    part_der_list_j = Jacobian_matrix_j(d, :);
                    
                    uncert_prop_matrix_ij = (sigma_list_i' * sigma_list_j) .* (part_der_list_i' * part_der_list_j) .* correlation_coeff_matrix_ij;

                    % The effect of each dimension of the variables is added
                    uncert_prop_ij              = sum(uncert_prop_matrix_ij, 'all');
                    uncert_prop_matrix(i, j, d) = uncert_prop_ij;
                    uncert_prop_matrix(j, i, d) = uncert_prop_ij;
                end
            end
        end
        
        % The total propagated uncertainty
        cylinder_centre_variance    = squeeze(sum(uncert_prop_matrix, [1, 2]))';        % The matrix is summed from n*n*3 to 1*1*3 and then transformed to 1*3
        cylinder_centre_sigma       = sqrt(abs(cylinder_centre_variance));              % Negative variance cannot exist in reality, however (minor) negative variance is possible when propagating uncertainty using the method of moments

        % The confidence interval
        cylinder_centre_confidence_interval = cylinder_centre_mu + [1; -1] * cylinder_centre_sigma;
        
        Cylinder_Geometry.Cylinder_centre = struct('mu', cylinder_centre_mu, 'sigma', cylinder_centre_sigma, 'confidence_interval', cylinder_centre_confidence_interval);

        % Finally, the height was defined to be equal from the centre
        height_top_mu = length_mu / 2;
        height_bot_mu = -length_mu / 2;
        
        height_top_confidence_interval = height_top_mu + [1; -1] * height_top_sigma;
        height_bot_confidence_interval = height_bot_mu + [1; -1] * height_bot_sigma;
        
        Cylinder_Geometry.Cylinder_height_top = struct('mu', height_top_mu, 'sigma', height_top_sigma, 'confidence_interval', height_top_confidence_interval);
        Cylinder_Geometry.Cylinder_height_bot = struct('mu', height_bot_mu, 'sigma', height_bot_sigma, 'confidence_interval', height_bot_confidence_interval);

    %% Progress function %%
        function ProgressUpdate(~)
            % Ensures that at most every percent is printed
            tick = tick + 1;    

            progress_last   = (tick - 1) / N * 100;
            progress        = tick / N * 100;

            if floor(progress) - floor(progress_last) >= 1
                fprintf('   Monte Carlo cylinder fitting progress: %i cylinders \n', tick);
            end            
        end
    
end
