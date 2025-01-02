% This script fits a circle to the given data by using the least-squares approach
% It requires x and y coordinates to be given, and spits out the centre of the circle (circle_x, circle_y) as well as its radius (circle_r)

% Three methods can be used: algebraic least-squares (LS), geometric least-squares (LS-Geom) and Pratt's least-squares (LS-Pratt)
% They can be weighted by uncertainty if Weighting is true. This requires the locations of the scanners, laser beam divergence and range uncertainty of the scanner itself
% Alternatively, they can be simply weighted by distance too
% Note however that weighting is not currently implemented for geometric least-squares

function [circle_x, circle_y, circle_r] = Least_Squares_Circle_Fitting(x_cell, y_cell, Method, Weighting, Scanner_loc_cell, beam_divergence, sigma_range_device, range_bias, max_incidence_angle, Point_Duplication)
    
    %% Inputs %%
        % Threshold for convergence
        min_step_threshold  = 1e-6;         % [%] If the step size is below this value, it has converged
        max_iterations      = 1e2;          % [-] If the maximum number of iterations is exceeded the last values will be used and a message is printed

        Weighting_type      = 'Distance';   % [Distance, Uncertainty] Either can be used to weigh the data points
        
    %% The range bias is compensated for %%
        if range_bias ~= 0
            number_scanners = length(Scanner_loc_cell);

            for s = 1 : number_scanners
                % The coordinates are moved closer or farther away from the scanner
                scanner_loc     = Scanner_loc_cell{s};

                x_list          = x_cell{s};
                y_list          = y_cell{s};

                % The vectors from the scanner to the points
                vector_x_list   = x_list - scanner_loc(1);
                vector_y_list   = y_list - scanner_loc(2);
                
                % The vectors are normalised
                vector_n_list           = [vector_x_list, vector_y_list] ./ sqrt(vector_x_list.^2 + vector_y_list.^2);
                
                % The bias is compensated for
                bias_compensation_list  = range_bias * vector_n_list;
                
                x_list_comp = x_list - bias_compensation_list(:, 1);
                y_list_comp = y_list - bias_compensation_list(:, 2);
                
                x_cell{s}   = x_list_comp;
                y_cell{s}   = y_list_comp;
            end
        end
        
    %% Points from all scanners are combined %%        
        x_list = vertcat(x_cell{:});
        y_list = vertcat(y_cell{:});
        
        % A column vector shape is ensured
        n = length(x_list);

        x_list = reshape(x_list, [n, 1]);
        y_list = reshape(y_list, [n, 1]);
        
    %% The weights are computed beforehand if they are distance-based %%
        if strcmp(Weighting_type, 'Distance')
            % The average distance of each scanner to their points is computed and used for the weight of its points
            number_scanners         = length(Scanner_loc_cell);

            weights_cell            = cell(1, number_scanners);

            for s = 1 : number_scanners
                % The location of the scanner
                scanner_x       = Scanner_loc_cell{s}(1);
                scanner_y       = Scanner_loc_cell{s}(2);

                % This scanner's point cloud
                x_points_list   = x_cell{s};
                y_points_list   = y_cell{s};

                number_points   = length(x_points_list);

                % The distance between the scanner and each point in its point cloud
                scanner_distance_list   = sqrt((scanner_x - x_points_list).^2 + (scanner_y - y_points_list).^2);

                % The average distance is used
                scanner_distance        = mean(scanner_distance_list);

                % The weight of these points is the inverse of their distance and equal for all points
                weight          = 1 / scanner_distance;
                weights_list    = weight * ones(number_points, 1);

                weights_cell{s} = weights_list;
            end

            w_list = vertcat(weights_cell{:});
        end
    
    %% Algebraic least squares %%
    if strcmp(Method, 'LS')
        %--% An unweighted fit is always performed first %--%
        % The least squares fit matrices (A*x = B)
        A = [x_list, y_list, ones(n, 1)];
        B = -(x_list.^2 + y_list.^2);        
        
        % The least squares solution
        circle_parameters = A \ B;

        circle_x = -circle_parameters(1) / 2;
        circle_y = -circle_parameters(2) / 2;
        circle_r = sqrt(circle_x^2 + circle_y^2 - circle_parameters(3));

        %--% Weighted fits %--%
        if  Weighting == true
            if strcmp(Weighting_type, 'Distance')                
                W = [w_list, w_list, w_list];

                % The matrices are changed accordingly
                C = A' * (W .* A);
                D = A' * (W .* B);

                % The least squares solution
                circle_parameters_w = C \ D;

                circle_x = -circle_parameters_w(1) / 2;
                circle_y = -circle_parameters_w(2) / 2;
                circle_r = sqrt(circle_x^2 + circle_y^2 - circle_parameters_w(3));
            
            elseif strcmp(Weighting_type, 'Uncertainty')
                iter = 0;

                circle_x_old = circle_x;
                circle_y_old = circle_y;
                circle_r_old = circle_r;

                while iter <= max_iterations
                    iter = iter + 1;

                    % The uncertainty of the entire point cloud
                    [sigma_radial_cell, sigma_range_cell, ~, ~] = Circular_Object_Uncertainty(circle_x_old, circle_y_old, circle_r_old, Scanner_loc_cell, x_cell, y_cell, beam_divergence, sigma_range_device, max_incidence_angle, Point_Duplication);

                    sigma_radial_total  = vertcat(sigma_radial_cell{:});
                    sigma_range_total   = vertcat(sigma_range_cell{:});

                    % The weights are computed using the uncertainty
                    w_list = sigma_radial_total .* sigma_range_total;

                    W = [w_list, w_list, w_list];

                    % The matrices are changed accordingly
                    C = A' * (W .* A);
                    D = A' * (W .* B);

                    % The least squares solution
                    circle_parameters_w = C \ D;

                    circle_x = -circle_parameters_w(1) / 2;
                    circle_y = -circle_parameters_w(2) / 2;
                    circle_r = sqrt(circle_x^2 + circle_y^2 - circle_parameters_w(3));

                    % Check for convergence
                    [Convergence, circle_x_old, circle_y_old, circle_r_old] = Convergence_Check(min_step_threshold, [circle_x_old, circle_y_old, circle_r_old], [circle_x, circle_y, circle_r]);

                    if Convergence == true
                        break
                    end
                end
            end
        end
        
    %% Geometrical least squares %%
    elseif strcmp(Method, 'LS-Geom')
        %--% An unweighted fit is performed for the initial guess %--%
        % The least squares fit matrices (A*x = B)
        A = [x_list, y_list, ones(n, 1)];
        B = -(x_list.^2 + y_list.^2);        
        
        % The least squares solution
        circle_parameters = A \ B;

        circle_x = -circle_parameters(1) / 2;
        circle_y = -circle_parameters(2) / 2;
        circle_r = sqrt(circle_x^2 + circle_y^2 - circle_parameters(3));        
        
        % Gauss-Newton method to solve the non-linear least-squares problem
        convergence = false;
        iter        = 0;
        
        while convergence == false && iter < max_iterations
            iter = iter + 1;
            
            % The residuals
            R = sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2) - circle_r;

            % The derivatives
            dRdx_g = (circle_x - x_list) ./ sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2);
            dRdy_g = (circle_y - y_list) ./ sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2);
            dRdr_g = -ones(n, 1);

            % The Jacobian
            J = [dRdx_g, dRdy_g, dRdr_g];

            % The next geometry parameters
            Beta        = [circle_x; circle_y; circle_r];
            Beta_new    = Beta - (J' * J)\(J' * R);

            circle_x    = Beta_new(1);
            circle_y    = Beta_new(2);
            circle_r    = Beta_new(3);

            % Check for convergence of the squared sum of residuals
            sum_residuals = sum(R.^2);
            
            if iter == 1
                sum_residuals_init = sum_residuals;
            end

            R_new = sqrt((x_list - circle_x).^2 + (y_list - circle_y).^2) - circle_r;        
            sum_residuals_new = sum(R_new.^2);

            change = (sum_residuals_new - sum_residuals)/sum_residuals;
            
            % Alternatively, if the relative step size is low it is also said to have converged
            relative_step_size = abs(sum_residuals_new - sum_residuals) / sum_residuals_init;
                        
            if abs(change) < min_step_threshold
                convergence = true;
            elseif relative_step_size < min_step_threshold
                convergence = true;
            end            
        end
                        
    %% Pratt's least-squares approach for circle fitting %%
    elseif strcmp(Method, 'LS-Pratt')
        if Weighting == false
            [circle_x, circle_y, circle_r] = Pratt_Circle_Fit(x_list, y_list, min_step_threshold, max_iterations, Weighting, []);      
        
        else
            if strcmp(Weighting_type, 'Distance')
                [circle_x, circle_y, circle_r] = Pratt_Circle_Fit(x_list, y_list, min_step_threshold, max_iterations, Weighting, w_list);      

            elseif strcmp(Weighting_type, 'Uncertainty')
                %--% An unweighted fit is always performed first such that the uncertainty (and thus weights) can be computed %--%
                [circle_x, circle_y, circle_r] = Pratt_Circle_Fit(x_list, y_list, min_step_threshold, max_iterations, false, []);
        
                %--% Weighted fits are performed iteratively, as the uncertainty changes each time %--%
                iter = 0;

                circle_x_old = circle_x;
                circle_y_old = circle_y;
                circle_r_old = circle_r;

                while iter <= max_iterations
                    iter = iter + 1;

                    % The uncertainty of the entire point cloud
                    [sigma_radial_cell, sigma_range_cell, ~, ~] = Circular_Object_Uncertainty(circle_x_old, circle_y_old, circle_r_old, Scanner_loc_cell, x_cell, y_cell, beam_divergence, sigma_range_device, max_incidence_angle, Point_Duplication);

                    sigma_radial_total  = vertcat(sigma_radial_cell{:});
                    sigma_range_total   = vertcat(sigma_range_cell{:});

                    % The weights are computed based on the uncertainty
                    w_list = sigma_radial_total .* sigma_range_total;

                    [circle_x, circle_y, circle_r] = Pratt_Circle_Fit(x_list, y_list, min_step_threshold, max_iterations, Weighting, w_list);      

                    % Check for convergence
                    [Convergence, circle_x_old, circle_y_old, circle_r_old] = Convergence_Check(min_step_threshold, [circle_x_old, circle_y_old, circle_r_old], [circle_x, circle_y, circle_r]);

                    if Convergence == true
                        break
                    end
                end
            end
        end
    end

%% Convergence check function %%
function [Convergence, circle_x_old, circle_y_old, circle_r_old] = Convergence_Check(Step_threshold, old_geometry, new_geometry)
    % Check for convergence
    step_sizes      = (new_geometry - old_geometry) ./ old_geometry * 100;

    % The largest step size is used
    step = max(abs(step_sizes));

    if abs(step) < Step_threshold
        Convergence = true;

        circle_x_old = old_geometry(1);
        circle_y_old = old_geometry(2);
        circle_r_old = old_geometry(3);

    % If convergence is not reached, the old geometry is updated to the new one
    else
        Convergence = false;

        circle_x_old = new_geometry(1);
        circle_y_old = new_geometry(2);
        circle_r_old = new_geometry(3);
    end
end
end