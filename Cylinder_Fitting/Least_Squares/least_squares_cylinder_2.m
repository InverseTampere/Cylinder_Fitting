% This is an updated version of the script created by Pasi Raumonen for TreeQSM
% It performs least-squares cylinder using the Gauss-Newton method

% ---------------------------------------------------------------------
% Input    
    % point_cloud_matrix                Point cloud
    % Initial_Cylinder_Geometry         Initial estimates of the cylinder parameters
    % weight_list                       (Optional) Weights of the points for fitting
    % point_cloud_subset                (Optional) Subset of "P" where the cylinder is intended
 
% Output  
    % Cylinder_Geometry                 Structure array containing the following fields:
    %   radius                          Radius of the cylinder
    %   length                          Length of the cylinder
    %   start                           Point on the axis at the bottom of the cylinder (1 x 3)
    %   axis                            Axis direction of the cylinder (1 x 3) 
    %   mad                             Mean absolute distance between points and cylinder surface
    %   SurfCov                         Relative cover of the cylinder's surface by the points 
    %   dist                            Radial distances from the points to the cylinder (m x 1) 
    %   conv                            If conv = 1, the algorithm has converged 
    %   rel                             If rel = 1, the algorithm has reliable answer in terms of
    %                                   matrix inversion with a good enough condition number
% ---------------------------------------------------------------------

%%%%%% CHANGES %%%%%%
% Made parameter names more descriptive
% Changed convergence to be based on a relative difference, with the threshold given at the start of the file
% Introduced a maximum number of iterations warning
% Outputs are in double format now
% The mean absolute distance is now taken over all weighted points, not just the one with the highest weight
% The geometry variables are normalised when the Jacobian is computed
% The history of the variables is now saved to halve the number of computations and create diagnostic plots
% Changed from nargin to ~isempty to make it a bit more flexible
% Poor reliability often made the geometry update explode, so changed how the script deals with that

function Cylinder_Geometry = least_squares_cylinder_2(point_cloud_matrix, Initial_Cylinder_Geometry, weight_list, point_cloud_subset)

    %% Initial cylinder geometry %%
        cyl_dir_init    = Initial_Cylinder_Geometry.axis;
        cyl_start_init  = Initial_Cylinder_Geometry.start;
        cyl_radius_init = Initial_Cylinder_Geometry.radius;

    %% Manual inputs %%
        resolution              = 0.03;     % [m] "Resolution level" for computing the surface coverage
        num_max_iter            = 1e2;      % [-] Maximum number of Gauss-Newton iterations
        convergence_threshold   = 1e-3;     % [%] If the difference in the sum of squared distances between updates is below this value, the GN method is said to have converged

        Normalisation           = false;    % [true, false] If the variables are normalised or not
        bounds_margin           = 0.25;     % [-] Factor used for normalisation

        History_Plot            = false;    % [true, false] Plots the history of the variables during the GN algorithm

    %% Point cloud transformation %%
        % Transform the data to close to standard position via a translation  
        % followed by a rotation
        init_rotation_matrix    = rotate_to_z_axis(cyl_dir_init);
        point_cloud_matrix_t    = (point_cloud_matrix - cyl_start_init) * init_rotation_matrix';
        
        % Indices for the geometry vector to ensure consistency
        [centre_ind, alpha_ind, beta_ind, radius_ind]   = deal(1:2, 3, 4, 5);     
        Geometry_Indices                                = struct('centre', centre_ind, 'alpha', alpha_ind, 'beta', beta_ind, 'radius', radius_ind);
        num_geom_parameters                             = radius_ind;

        % Resulting initial geometry vector has only a nonzero radius (centre is the point the cylinder axis intersects the x-y plane, alpha is rotation around x-axis, beta around y-axis)
        cylinder_geometry       = zeros(num_geom_parameters, 1);
        cylinder_geometry([centre_ind, alpha_ind, beta_ind, radius_ind]) = [0, 0, 0, 0, cyl_radius_init]; 

    %% Variable normalisation %%
        % The cylinder variable's are normalised s.t. they (and their Jacobians) are of the same order of magnitude
        if Normalisation == true
            proj_centre_bounds  = bounds_margin * cyl_radius_init * [1, 1];
            angle_bounds        = bounds_margin * pi * [1, 1];
            radius_LB           = cyl_radius_init / (1 + bounds_margin);        % Note that because the radius is always positive, the radius bounds are not exactly symmetrical
            radius_UB           = cyl_radius_init * (1 + bounds_margin);
    
            geometry_LB         = zeros(num_geom_parameters, 1);
            geometry_UB         = zeros(num_geom_parameters, 1);
            geometry_LB([centre_ind, alpha_ind, beta_ind, radius_ind]) = [-proj_centre_bounds, -angle_bounds, radius_LB]';
            geometry_UB([centre_ind, alpha_ind, beta_ind, radius_ind]) = [proj_centre_bounds, angle_bounds, radius_UB]';

        else
            geometry_LB = zeros(num_geom_parameters, 1);        % Values of 0 and 1 mean that no change takes place
            geometry_UB = ones(num_geom_parameters, 1);
        end

    %% Gauss-Newton algorithm %%
        % The cylinder variables are found iteratively
        iter            = 0; 
        Convergence     = false;            % Did the iterations converge
        Reliability     = true;             % Did the problem become uninvertible

        % The history of the variables is saved for later reference
        SSD_list        = zeros(num_max_iter + 1, 1);                       % Sum of squared distances
        geometry_matrix = zeros(num_max_iter + 1, num_geom_parameters);     % Geometry variables

        % Values for the initial geometry
        if isempty(weight_list)
            [distance_list, SSD, Jacobian] = func_grad_cylinder_2(cylinder_geometry, Geometry_Indices, point_cloud_matrix_t);
        else
            [distance_list, SSD, Jacobian] = func_grad_cylinder_2(cylinder_geometry, Geometry_Indices, point_cloud_matrix_t, weight_list);
        end

        SSD_list(1)             = SSD;
        geometry_matrix(1, :)   = cylinder_geometry;

        while iter <= num_max_iter && ~Convergence               
            %--% Update step %--%
            % Update the iteration counter
            iter = iter + 1;

            % The normalised Jacobian
            Jacobian_n = Jacobian .* (geometry_UB - geometry_LB)';
  
            % solve for the system of equations: par(i+1) = par(i) - (J'J)^(-1)*J'distance_list(par(i))
            A = Jacobian_n'*Jacobian_n;
            b = Jacobian_n'*distance_list;

            % Invertability check
            condition_number_A = rcond(-A);

            if condition_number_A < 1e4 * eps        
                % In this case the problem is ill-suited to inversion and the scheme stops - as otherwise the update step likely explodes
                geometry_update_n   = zeros(num_geom_parameters, 1);                % Note that this causes the convergence to be true
                Reliability         = false;
            else
                % The system of equations produces the normalised geometry update
                % Warnings are turned off for the condition number, as it is evaluated differently
                warning off
                geometry_update_n   = -A\b;     
                warning on      
            end

            % If the update is well out of normalised bounds, the matrix operation is assumed to have failed
            % The update step is set to 0, forcing convergence
            if max(abs(geometry_update_n)) > 1e1            % Note that they should be normalised between 0 and 1
                geometry_update_n = zeros(num_geom_parameters, 1);
            end
            
            % The updated normalised geometry vector
            cylinder_geometry_n = (cylinder_geometry - geometry_LB) ./ (geometry_UB - geometry_LB);
            cylinder_geometry_n = cylinder_geometry_n + geometry_update_n;

            % A check is used to ensure the updated geometry is within the bounds
            cylinder_geometry_n = max(0, cylinder_geometry_n);
            cylinder_geometry_n = min(1, cylinder_geometry_n);

            cylinder_geometry = cylinder_geometry_n .* (geometry_UB - geometry_LB) + geometry_LB;

            geometry_matrix(iter + 1, :) = cylinder_geometry;                                       % And saved for future reference          

            %--% SSD and Jacobian for the updated geometry %--%
            if isempty(weight_list)
                [distance_list, SSD, Jacobian] = func_grad_cylinder_2(cylinder_geometry, Geometry_Indices, point_cloud_matrix_t);
            else
                [distance_list, SSD, Jacobian] = func_grad_cylinder_2(cylinder_geometry, Geometry_Indices, point_cloud_matrix_t, weight_list);
            end
            
            SSD_list(iter + 1) = SSD;

            %--% Convergence check %--%
            % Relative change
            SSD_prev            = SSD_list(iter);                           % Squared sum of the distances of the geometry before the update
            rel_distance_change = abs((SSD - SSD_prev) / SSD_prev) * 100;   % Relative change in percent

            if rel_distance_change < convergence_threshold
                Convergence = true;
            end
        end

        % A warning is displayed if the maximum number of iterations was reached
        if iter == num_max_iter
            warning('The maximum number of iterations during least-squares fitting was reached.');

        % If it wasn't, the history is trimmed
        else
            SSD_list        = SSD_list(1 : iter + 1);
            geometry_matrix = geometry_matrix(1 : iter + 1, :);
        end

        % If the reliability was false, convergence is also made false now
        if Reliability == false
            Convergence = false;
        end

        % Fitted geometry
        [proj_centre, alpha, beta, cyl_radius] = deal(cylinder_geometry(centre_ind), cylinder_geometry(alpha_ind), cylinder_geometry(beta_ind), cylinder_geometry(radius_ind));

        % History of the Gauss-Newton algorithm
        if History_Plot == true
            % Number of rows and columns for the subplot
            number_rows     = 2;
            number_columns  = ceil((num_geom_parameters + 1) / number_rows);

            % Colour map
            colour_map = cbrewer('qual', 'Set1', num_geom_parameters + 1);

            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on

            %--% Geometry variables %--%
            geom_labels     = fieldnames(Geometry_Indices);
            dim_labels      = {'x', 'y'};
            num_geom_labels = length(geom_labels);

            for g = 1 : num_geom_labels
                geom_label  = geom_labels{g};
                geom_ind    = Geometry_Indices.(geom_label);
                num_dim     = length(geom_ind);
                geom_colour = colour_map(g, :);

                for d = 1 : num_dim
                    % This geometry's data
                    ind         = geom_ind(d);
                    geom_data   = geometry_matrix(:, ind);
                    
                    if num_dim > 1
                        geom_label = geom_labels{g};
                        dim_label  = dim_labels{d};
                        geom_label = sprintf('%s_%s', geom_label, dim_label);
                    end

                    % Subplot
                    subplot(number_rows, number_columns, ind)
                    hold on
                    grid on

                    plot(geom_data, 'LineWidth', 2, 'Color', geom_colour);

                    % Axes
                    xlabel('Iteration [-]');
                    ylabel(sprintf('%s [m]', geom_label));
                    xlim([1, iter]);
                    
                    % Font size
                    set(gca, 'FontSize', 15);
                    set(gca, 'LineWidth', 2);
        
                    hold off
                end
            end

            %--% Distance %--%
            distance_colour = colour_map(num_geom_parameters + 1, :);

            subplot(number_rows, number_columns, num_geom_parameters + 1)
            hold on
            grid on

            plot(SSD_list, 'LineWidth', 2, 'Color', distance_colour);
            
            % Axes
            xlabel('Iteration [-]');
            ylabel('SSD [m]');
            xlim([1, iter]);

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            % Pause message
            disp('The least-squares cylinder fitting script will finish and the plot will close when a key is pressed');
            pause();

            close(1);             
        end

    %% Compute the cylinder parameters and other outputs %%
        % Cylinder axis
        rotation_matrix = form_rotation_matrices_2(alpha, beta);
        z_axis          = [0, 0, 1];                                                % Note that the original rotation matrix was defined to have the cylinder axis as the z-axis
        cyl_axis        = (init_rotation_matrix' * rotation_matrix' * z_axis')';

        % If the point cloud subset is given and is sufficiently large, it is used for the start, length, mad, and SurfCov
        if ~isempty(point_cloud_subset)
            num_subset_points = size(point_cloud_subset, 1);

            if num_subset_points > 5
                point_cloud_matrix = point_cloud_subset;
            end
        end

        % The length
        height_list = point_cloud_matrix * cyl_axis';               % heights along the axis
        height_min  = min(height_list);
        height_max  = max(height_list);

        cyl_length  = height_max - height_min;

        % Cylinder start point
        proj_centre_3D  = [proj_centre; 0]';
        axis_point      = (init_rotation_matrix' * proj_centre_3D')' + cyl_start_init;         % Point on axis
        height_point    = cyl_axis * axis_point';
        cyl_start       = axis_point - (height_point - height_min) * cyl_axis;              % Axis point at the cylinder's bottom    

        % Mean absolute distance for the points along the cylinder length
        MAD = mean(abs(distance_list));           

        % Compute surface coverage, minimum 3*8 grid
        if ~any(isnan(cyl_axis)) && ~any(isnan(axis_point)) && Convergence
            nl = max(3, ceil(cyl_length/resolution));
            ns = ceil(2*pi * cyl_radius/resolution);
            ns = min(36, max(ns, 8));
    
            SurfCov = surface_coverage(point_cloud_matrix, cyl_axis, axis_point, nl, ns, 0.8*cyl_radius);
        else
            SurfCov = 0;
        end

        % Outputs in structure format
        Cylinder_Geometry = struct('radius', cyl_radius, 'length', cyl_length, 'start', cyl_start, 'axis', cyl_axis, ... 
                                   'mad', MAD, 'SurfCov', SurfCov, 'dist', distance_list, 'conv', Convergence, 'rel', Reliability);

end
