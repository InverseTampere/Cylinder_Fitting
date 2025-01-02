% The Mahalanobis distance (NOT the expected distance!) is computed using the line approximation between the given distributions' mu and cylinder surface

function [avg_Mahal_distance, Mahal_distance_list] = Mahalanobis_Distance_Cylinder_Line_Approx(cylinder_radius, cylinder_centre, cylinder_direction, Point_Cloud_Distributions, Scanner_loc_cell, distance_moment, Plot)

    %% Point cloud projection onto cross-section %%
        [mu_matrix, sigmae_matrix]  = Distribution_to_Cylinder_Cross_Section_Projection(cylinder_centre, cylinder_direction, Point_Cloud_Distributions, Scanner_loc_cell);
        [number_points, ~]          = size(mu_matrix);

    %% Mahalanobis distance using the line approximation %%
        % The angle between the y-axis and vector form mu towards the origin (in the projected frame)
        beta_list       = atan(abs(mu_matrix(:, 1) ./ mu_matrix(:, 2)));        % [0, pi/2]
        oblique_bool    = beta_list > pi/2 - 1e-6;                              % A different approach is needed for oblique angles. A small margin is used for numerical rounding

        % The distance between the projected mu and the circle centre
        norm_mu_list    = sqrt(sum(mu_matrix.^2, 2));
        rho_list        = abs(norm_mu_list - cylinder_radius);
        s_list          = rho_list ./ cos(beta_list);                           % Infinite when beta = pi/2

        % The Mahalanobis distance
        Mahal_distance_list                 = s_list ./ sqrt(sigmae_matrix(:, 1).^2.*tan(beta_list).^2 + sigmae_matrix(:, 2).^2);
        Mahal_distance_list(oblique_bool)   = rho_list(oblique_bool) ./ sigmae_matrix(oblique_bool, 1);

        Mahal_distance_list = Mahal_distance_list.^distance_moment;         % The moment is taken
        avg_Mahal_distance  = mean(Mahal_distance_list);

    %% Plot %%
        if Plot == true
            % Distribution properties
            distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
            distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
            distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;

            % The number of coordinates used for each distribution and the cylinder
            number_coord = 1e2;     

            % Estimate of the length for the plot
            mu_0_matrix         = vertcat(distribution_mu_cell{:}) - cylinder_centre;
            mu_0_norm_list      = sqrt(sum(mu_0_matrix.^2, 2));
            cylinder_length     = 2*max(mu_0_norm_list);

            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
            
            % Cylinder
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Cylinder');

            % Distributions and their distances           
            for d = 1 : number_points
                % Distribution
                sigma_x     = sigmae_matrix(d, 1);          % Used to scale the Mahalanobis distance
                sigma_y     = sigmae_matrix(d, 2);

                mu          = distribution_mu_cell{d};
                sigmae      = distribution_sigmae_cell{d};
                distr_axes  = distribution_axes_cell{d};

                [distr_coord_matrix, number_coord_distr] = Ellipsoid_Coordinate_Generator(mu, sigmae, distr_axes, number_coord);

                x_ellipsoid = reshape(distr_coord_matrix(:, 1), sqrt(number_coord_distr) * [1, 1]);
                y_ellipsoid = reshape(distr_coord_matrix(:, 2), sqrt(number_coord_distr) * [1, 1]);
                z_ellipsoid = reshape(distr_coord_matrix(:, 3), sqrt(number_coord_distr) * [1, 1]);
            
                surf_distr  = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.10, 'DisplayName', '1 \sigma');
                sc_mu       = scatter3(mu(1), mu(2), mu(3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', '\mu');
                
                % Mahalanobis distance, adjusted to be in metres
                Mahal_dist  = Mahal_distance_list(d)^(1/distance_moment);

                [proj_mu, ~, ~] = Point_to_Vector_Projection(mu, cylinder_direction, cylinder_centre);

                vector          = mu - proj_mu;

                vector_sign     = sign(cylinder_radius - norm(vector));     % Positive if it points outward
                vector          = sqrt(sigma_x * sigma_y) * Mahal_dist * vector_sign * vector / norm(vector);
                
                if distance_moment == 1
                    vec_name = sprintf('M * sqrt(%s %s)', '\sigma_x', '\sigma_y');
                elseif distance_moment == 2
                    vec_name = sprintf('sqrt(M^2 * %s %s)', '\sigma_x', '\sigma_y');
                end
                
                pl_vec = plot3(mu(1) + [0, vector(1)], mu(2) + [0, vector(2)], mu(3) + [0, vector(3)], 'color', 'k', 'lineWidth', 2, 'DisplayName', vec_name);
                
                if d > 1
                    surf_distr.HandleVisibility = 'Off';
                    sc_mu.HandleVisibility      = 'Off';
                    pl_vec.HandleVisibility     = 'Off';
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
            legend('show', 'location', 'eastoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off
            
            fprintf('The mean M^%g = %.3g for the cylinder.', distance_moment, avg_Mahal_distance);
            disp('The script will finish and figure will close upon a key-press.');
            pause();
            
            close(1);
        end
end