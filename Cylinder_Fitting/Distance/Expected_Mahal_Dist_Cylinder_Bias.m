% The expected Mahalanobis distance results in positive radial bias
% Here it is quantified as the delta between the optimum and norm for each point respectively

function delta_radius_list = Expected_Mahal_Dist_Cylinder_Bias(circle_centre, Projected_Distributions, number_samples, alpha, Diagnostics)
    
    %% Point cloud distributions %%
        distr_mu_cell       = Projected_Distributions.Plane.Projection.mu;
        distr_sigmae_cell   = Projected_Distributions.Plane.Projection.sigmae;
        distr_axes_cell     = Projected_Distributions.Plane.Projection.distr_axes;

    %% Manual inputs %%
        % Golden-section search
        convergence_threshold   = 1e-3;         % [-] Expressed in the normalised expected Mahalanobis distance
        max_iterations          = 020;          % [-]
        GSS_Print               = false;        % [true, false] Shows possible warning messages
        GSS_Diagnostics         = false;        % [true, false] Presents results for diagnosis
                
    %% Quantifying the radial bias %%
        % The radial bias for each distribution
        Radius_Bias_fun     = @(mu, sigmae, distr_axes) Radius_Bias(mu, sigmae, distr_axes);
        delta_radius_list   = cellfun(Radius_Bias_fun, distr_mu_cell, distr_sigmae_cell, distr_axes_cell);

        % A local function is used over each of the distributions
        function delta_radius = Radius_Bias(mu, sigmae, distr_axes)
            % Structure containing the data for this distribution
            Projected_Distribution.Plane.Projection = struct('mu', {{mu}}, 'sigmae', {{sigmae}}, 'distr_axes', {{distr_axes}});
            
            % The expected distance at the norm of mu for normalisation
            mu_norm = norm(mu);
            expected_Mahal_dist_norm = Expected_Mahalanobis_Distance_Cylinder_Numerical(circle_centre, mu_norm, [], Projected_Distribution, number_samples, alpha);
    
            % Function handle for the normalised expected Mahalanobis distance as a function of radius
            Norm_Mahal_Distance_Circle_fun = @(radius) Expected_Mahalanobis_Distance_Cylinder_Numerical(circle_centre, radius, [], Projected_Distribution, number_samples, alpha) / expected_Mahal_dist_norm;

            % The radius for which it is minimal is searched for within the norm and the norm plus/minus one standard deviation
            % The true solution is always greater than the norm of mu, however the minus deviation improves convergence of GSS
            search_lower_bound      = max(mu_norm * 1e-2, mu_norm - max(sigmae));       % Note that the minimum radius must be greater than zero
            search_upper_bound      = mu_norm + max(sigmae);
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
                distr_coord_matrix  = Ellipse_Coordinate_Generator(mu, distr_axes, sigmae, number_samples);
                plot(distr_coord_matrix(:, 1), distr_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', sprintf('Distribution, 1 %s', '\sigma'));
                
                % Circle with norm of mu
                [x_circle_list, y_circle_list] = Circle_Coordinates(mu_norm, 0, 0, number_samples);
                plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Circle with r = ||\mu||');                
                
                % Optimal circle
                [x_circle_list, y_circle_list] = Circle_Coordinates(optimal_radius, 0, 0, number_samples);
                plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'b', 'LineStyle', ':', 'DisplayName', 'Circle with optimal radius');                
                
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