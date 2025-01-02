% The geometry's convexity affects the expected distance resulting from the 1D scanning uncertainty and 2D power uncertainty to be positively biased
% The resulting bias to a greater radius is quantified here, both separately and combined

function [prop_bias_matrix, delta_radius_list, scanning_radius_bias_list, power_radius_bias_list] = Expected_Distance_Cylinder_Bias(distance_moment, cylinder_radius, Projected_Distributions_cell, alpha, number_samples, Plot, Diagnostics)

    %% Manual inputs %%    
        % Sampling search for the convexity bias
        sigma_ratio                 = 0.50;     % [-] The balanced radius is searched for within ||mu|| +/- ratio*distribution extent
        convergence_threshold_rel   = 1e-2;     % [-] Threshold in expected distance relative to the uncertainty
        max_iterations              = 020;      % [-]
        Print_GSS                   = false;    % [true, false]
        Diagnostics_GSS             = false;    % [true, false]

    %% Scanning bias %%
        % Samples given the confidence interval and unit uncertainty
        m_STD_vector    = sqrt(chi2inv(1 - alpha, 1));
        unit_p_list     = m_STD_vector * linspace(-1, 1, number_samples);
        unit_dp         = 2 * m_STD_vector / (number_samples + 1);
        
        % Their density values
        unit_f_list     = 1/sqrt(2*pi) * exp(-1/2 * unit_p_list.^2);
        
        % The scanning bias for each projected distribution
        Scanning_Bias_fun           = @(Projected_Distributions) Scanning_Bias(Projected_Distributions); 
        scanning_radius_bias_list   = cellfun(Scanning_Bias_fun, Projected_Distributions_cell);
        
        function scanning_bias = Scanning_Bias(Projected_Distributions)
            % Global variables are assigned for variables that are needed for the minimum search function
            global step_size density_list norm_list
            
            % Projected distribution properties in the propagation direction
            mu          = Projected_Distributions.Plane.Projection.mu{1};
            sigma_p     = Projected_Distributions.Plane.Projection.sigmae{1}(2);
            prop_axis   = Projected_Distributions.Plane.Projection.distr_axes{1}(2, :);     
            
            % The samples, their density values and distance between samples get scaled accordingly
            p_list          = unit_p_list * sigma_p;
            density_list    = unit_f_list / sigma_p;
            step_size       = unit_dp * sigma_p;
            
            % The incidence angle is the acosine of the absolute dot product between the propagation axis and normal vector
            mu_norm         = norm(mu);
            normal_vector   = mu / mu_norm;
            incidence_angle = acos(abs(dot(normal_vector, prop_axis)));
            
            % The coordinates of mu are transformed to the beam-oriented but circle-centered frame
            mu_r = mu_norm * sin(incidence_angle);
            mu_p = mu_norm * -cos(incidence_angle);
            
            norm_list = sqrt(mu_r^2 + (mu_p + p_list).^2);
            
            % The radius for which the expected distance is minimal is found iteratively
            search_lower_bound  = mu_norm - sigma_p * sigma_ratio;
            search_upper_bound  = mu_norm + sigma_p * sigma_ratio;
            search_lower_bound  = max(1e-6, search_lower_bound);            % Note that the radius must be > 0       
            
            convergence_threshold       = convergence_threshold_rel * sigma_p^distance_moment;
            
            Expected_Distance_fun       = @(radius) Expected_Distance(radius);
            [minimum_ED_radius, ~, ~]   = Golden_Section_Search(Expected_Distance_fun, search_lower_bound, search_upper_bound, convergence_threshold, max_iterations, Print_GSS, Diagnostics_GSS);
            
            % The delta to the norm of mu
            scanning_bias = minimum_ED_radius - mu_norm;
            
            % Diagnostics plot to see whether the minimum was found correctly
            if Diagnostics == true
                search_radius_list      = linspace(search_lower_bound, search_upper_bound, number_samples);
                expected_distance_list  = arrayfun(Expected_Distance_fun, search_radius_list);
                
                figure(1)
                % Size and white background
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    

                hold on
                grid on
            
                plot(search_radius_list, expected_distance_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('E[D^%g(r)', distance_moment));
                plot([minimum_ED_radius, minimum_ED_radius], [0, 1.2*max(expected_distance_list)], 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Found minimum');
                plot([mu_norm, mu_norm], [0, 1.2*max(expected_distance_list)], 'LineWidth', 2, 'color', 'r', 'DisplayName', '||\mu||');
                
                 % Axes
                xlabel('r [m]')
                ylabel(sprintf('E[D^%g] [m^%g]', distance_moment, distance_moment))
                ylim([0, 1.2*max(expected_distance_list)]);
                
                % Legend
                legend('show', 'location', 'eastoutside');

                % Font size
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            
                % Pause
                disp('Diagnostics are shown for the radius/E[D] relation The figure will close and script will continue upon a key-press.');
                pause();

                close(1);
            end
        end
    
    %% Power bias %%            
        % This embedded function finds the delta between the expected geometry and that which minimises the expected distance for each distribution
        Power_Bias_fun                      = @(Projected_Distributions) Power_Bias(Projected_Distributions);
        [prop_bias_cell, radius_bias_cell]  = cellfun(Power_Bias_fun, Projected_Distributions_cell, 'UniformOutput', false); 
                  
        prop_bias_matrix        = vertcat(prop_bias_cell{:});
        power_radius_bias_list  = vertcat(radius_bias_cell{:});
    
        function [propagation_bias, radius_bias] = Power_Bias(Projected_Distributions)
            % Global variables are assigned for variables that are needed for the minimum search function
            global step_size density_list norm_list
            
            % Projected distribution properties
            mu      = Projected_Distributions.Plane.Original.mu{1};
            sigmae  = Projected_Distributions.Plane.Original.sigmae{1};
            axes    = Projected_Distributions.Plane.Original.distr_axes{1};
            
            % The incidence angle of the beam w.r.t. the circle
            mu_norm         = norm(mu);
            normal_vector   = mu / mu_norm;
            
            prop_axis       = axes(2, :);
            incidence_angle = acos(abs(dot(normal_vector, prop_axis)));
            
            % The expected value and covariance matrix are transformed to the beam-oriented but circle-centered frame
            mu_r    = mu_norm * sin(incidence_angle);
            mu_p    = mu_norm * -cos(incidence_angle);
            mu_rp   = [mu_r, mu_p];
            
            identity_matrix = eye(2);
            Sigma_rp        = sigmae .* identity_matrix;
            
            % Samples within the confidence interval of the distribution
            m_STD_ellipse   = sqrt(chi2inv(1 - alpha, 2));
            radii           = m_STD_ellipse * sigmae;
            [distr_samples_matrix, number_samples_adj] = Equal_Area_Elliptical_Sampler(mu_rp, radii, identity_matrix, number_samples);
        
            distr_area  = pi * prod(radii);
            step_size   = distr_area / number_samples_adj;
            
            % Their density and norm
            density_list    = mvnpdf(distr_samples_matrix, mu_rp, Sigma_rp);
            norm_list       = sqrt(sum(distr_samples_matrix.^2, 2));
            
            % The most optimal propagation coordinate for the circle is its expected value location along that axis
            propagation_bias    = mu_p * prop_axis;
            
            % The radius for which the expected distance is minimal is found iteratively
            sigma_r             = sigmae(1);
            search_lower_bound  = mu_norm - sigma_r * sigma_ratio;
            search_upper_bound  = mu_norm + sigma_r * sigma_ratio;
            search_lower_bound  = max(1e-6, search_lower_bound);            % Note that the radius must be > 0       
            
            convergence_threshold       = convergence_threshold_rel * sigma_r^distance_moment;
            
            Expected_Distance_fun       = @(radius) Expected_Distance(radius);
            [minimum_ED_radius, ~, ~]   = Golden_Section_Search(Expected_Distance_fun, search_lower_bound, search_upper_bound, convergence_threshold, max_iterations, Print_GSS, Diagnostics_GSS);
            
            % The delta to the radial coordinate of mu
            radius_bias = minimum_ED_radius - mu_r;
            
            % Diagnostics plot to see whether the minimum was found correctly
            if Diagnostics == true
                search_radius_list      = linspace(search_lower_bound, search_upper_bound, number_samples);
                expected_distance_list  = arrayfun(Expected_Distance_fun, search_radius_list);
                
                figure(1)
                % Size and white background
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    

                hold on
                grid on
            
                plot(search_radius_list, expected_distance_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('E[D^%g(r)', distance_moment));
                plot([minimum_ED_radius, minimum_ED_radius], [0, 1.2*max(expected_distance_list)], 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Found minimum');
                plot([mu_norm, mu_norm], [0, 1.2*max(expected_distance_list)], 'LineWidth', 2, 'color', 'r', 'DisplayName', '||\mu||');
                
                 % Axes
                xlabel('r [m]')
                ylabel(sprintf('E[D^%g] [m^%g]', distance_moment, distance_moment))
                ylim([0, 1.2*max(expected_distance_list)]);
                
                % Legend
                legend('show', 'location', 'eastoutside');

                % Font size
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            
                % Pause
                disp('Diagnostics are shown for the radius/E[D] relation The figure will close and script will continue upon a key-press.');
                pause();

                close(1);
            end
        end
  
    %% Total bias %%
        % The total bias is then the combined result
        delta_radius_list = scanning_radius_bias_list + power_radius_bias_list;
    
    %% Plot %%
        if Plot == true
            number_coord = 1e2;     % Number of coordinates used for the circle and each projected distribution
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
            
            % The circle
            [x_circle_list, y_circle_list] = Circle_Coordinates(cylinder_radius, 0, 0, number_coord);
            plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Circle');

            % The projected point cloud
            num_distributions = length(Projected_Distributions_cell);
            
            for i = 1 : num_distributions
                % Projected distribution
                Projected_Distributions = Projected_Distributions_cell{i};
                
                mu      = Projected_Distributions.Plane.Projection.mu{1};
                sigmae  = Projected_Distributions.Plane.Projection.sigmae{1};
                axes    = Projected_Distributions.Plane.Projection.distr_axes{1};
                
                distr_coord_matrix = Ellipse_Coordinate_Generator(mu, axes, sigmae, number_coord);
                
                pl_distr = plot(distr_coord_matrix(:, 1), distr_coord_matrix(:, 2), 'LineWidth', 1, 'color', 'r', 'DisplayName', sprintf('Point cloud, 1 %s', '\sigma'));
                sc_distr = scatter(mu(1), mu(2), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', sprintf('%s', '\mu'));
                             
                % Total radial bias and equivalent mu location
                delta_radius    = delta_radius_list(i);
                vector          = mu / norm(mu);
                
                sc_delta_distr  = scatter(mu(1) + delta_radius*vector(1), mu(2) + delta_radius*vector(2), 'filled', 'MarkerFaceColor', 'b', 'DisplayName', sprintf('%s r %s', '\Delta', '\mu'));
                pl_delta_rad    = plot(mu(1) + [0, delta_radius*vector(1)], mu(2) + [0, delta_radius*vector(2)], 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('%s r', '\Delta'));
                
                % Scanning bias
                scanning_bias   = scanning_radius_bias_list(i);                
                pl_scan_bias    = plot(mu(1) + [0, scanning_bias*vector(1)], mu(2) + [0, scanning_bias*vector(2)], 'LineWidth', 1, 'color', 'm', 'DisplayName', 'Scanning bias');
                                
                % Power bias
                radius_bias     = power_radius_bias_list(i); 
                prop_bias       = prop_bias_matrix(i, :);
                pl_power_bias   = plot(mu(1) + [0, radius_bias*vector(1) + prop_bias(1)], mu(2) + [0, radius_bias*vector(2) + prop_bias(2)], 'LineWidth', 1, 'color', 'c', 'DisplayName', 'Power bias');
                                                
                if i > 1
                    pl_distr.HandleVisibility       = 'Off';
                    sc_distr.HandleVisibility       = 'Off';
                    sc_delta_distr.HandleVisibility = 'Off';
                    pl_delta_rad.HandleVisibility   = 'Off';
                    pl_scan_bias.HandleVisibility   = 'Off';
                    pl_power_bias.HandleVisibility = 'Off';
                end
            end
            
            % Aspect ratio
            axis equal

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            
            % Legend
            legend('show', 'location', 'eastoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off
            
            % Pause
            disp('Expected distance bias quantification has finished, and the figure will close and script will continue upon a key-press.');
            pause();
            
            close(1);
        end
   
    %% Expected distance to the circle %%
        % Used to find the radius corresponding to the minimum expected distance
        function expected_distance = Expected_Distance(radius)
            % The distribution's density values and distances to the origin are retrieved
            global step_size density_list norm_list
            
            % The (expected) distance
            distance_list       = abs(norm_list - radius).^distance_moment;                         % Distance to the desired moment
            expected_distance   = sum(density_list .* distance_list) * step_size / (1 - alpha);     % Note that the finite confidence interval width is taken into account
        end
        
end