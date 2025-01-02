% The expected Mahalanobis distance results in positive radial bias
% Here it is quantified as the delta between the optimum and norm for each point respectively

function delta_radius_list = Expected_Mahal_Dist_Circle_Bias(circle_centre, Point_Cloud_Distributions, number_samples, Diagnostics)
    
    %% Point cloud distributions %%
        distr_mu_cell       = Point_Cloud_Distributions.distribution_mu_cell;
        distr_sigmae_cell   = Point_Cloud_Distributions.distribution_sigmae_cell;
        distr_axes_cell     = Point_Cloud_Distributions.distribution_axes_cell;

    %% Manual inputs %%
        % Golden-section search
        convergence_threshold   = 1e-2;         % [-] Expressed in the normalised expected Mahalanobis distance
        max_iterations          = 1e2;          % [-]
        GSS_Print               = true;        % [true, false] Shows possible warning messages
        GSS_Diagnostics         = true;        % [true, false] Presents results for diagnosis
                
    %% Translating and rotating the distributions %%
        % The point cloud is centered on the circle
        Centering_fun   = @(mu) mu - circle_centre;
        distr_mu_cell_c = cellfun(Centering_fun, distr_mu_cell, 'UniformOutput', false);
        
        % Rotated s.t. the distribution's axes are on the horizontal and vertical
        Rotation_fun    = @(mu, distr_axes) mu / distr_axes;
        distr_mu_cell_r = cellfun(Rotation_fun, distr_mu_cell_c, distr_axes_cell, 'UniformOutput', false);
        
    %% Quantifying the radial bias %%
        % The radial bias for each distribution
        Radius_Bias_fun     = @(mu, sigmae) Radius_Bias(mu, sigmae);
        delta_radius_list   = cellfun(Radius_Bias_fun, distr_mu_cell_r, distr_sigmae_cell);

        % The used scaling factor
        distr_mu_matrix_r   = vertcat(distr_mu_cell_r{:});
        distr_sigmae_matrix = vertcat(distr_sigmae_cell{:});

        mu_norm_list        = sqrt(sum(distr_mu_matrix_r.^2, 2));
        omega               = (max(mu_norm_list) - min(mu_norm_list)) / min(distr_sigmae_matrix, [], 'all');        % The amplitude of distance to the circle centre relative to the minimum uncertainty        
        
        % A local function is used over each of the distributions
        function delta_radius = Radius_Bias(mu, sigmae)
            % The expected distance at the norm of mu for normalisation
            mu_norm = norm(mu);
            num_dim = length(mu);
            origin  = zeros(1, num_dim);    % Note that the distributions have been centered on the circle
            expected_Mahal_dist_norm = Expected_Mahalanobis_Distance_Circle(origin, mu_norm, mu, sigmae, omega, number_samples);
    
            % Function handle for the normalised expected Mahalanobis distance as a function of radius
            Norm_Mahal_Distance_Circle_fun = @(radius) Expected_Mahalanobis_Distance_Circle(origin, radius, mu, sigmae, omega, number_samples) / expected_Mahal_dist_norm;
            
            % The radius for which it is minimal is searched for within the norm and the norm plus one standard deviation
            search_lower_bound      = mu_norm - 3*max(sigmae);
            search_upper_bound      = mu_norm + 3*max(sigmae);
            [optimal_radius, ~, ~]  = Golden_Section_Search(Norm_Mahal_Distance_Circle_fun, search_lower_bound, search_upper_bound, convergence_threshold, max_iterations, GSS_Print, GSS_Diagnostics);
            
            % The radial bias is then the delta between optimum and norm of mu
            delta_radius = optimal_radius - mu_norm;
            
            % Diagnostics plot
            if Diagnostics == true
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])     

                hold on
                grid on
                
                % The distribution
                distr_axes          = eye(num_dim);     % Note that the orientation w.r.t. the circle does not matter
                distr_coord_matrix  = Ellipse_Coordinate_Generator(mu, distr_axes, sigmae, number_samples);
                plot(distr_coord_matrix(:, 1), distr_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', sprintf('Distribution, 1 %s', '\sigma'));
                
                % Circle with norm of mu
                [x_circle_list, y_circle_list] = Circle_Coordinates(mu_norm, 0, 0, number_samples);
                plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Circle with r = ||\mu||');                
                
                % Optimal circle
                [x_circle_list, y_circle_list] = Circle_Coordinates(optimal_radius, 0, 0, number_samples);
                plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Circle with optimal radius');                
                
                % Axes
                axis equal
                xlabel('x [m]');
                ylabel('y [m]');

                % Legend
                legend('show', 'location', 'northoutside');

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off         

                % Pause
                disp('The line has been fitted. The figure will close and script will finish upon a key-press.');
                pause();

                close(1);
            end
        end

end