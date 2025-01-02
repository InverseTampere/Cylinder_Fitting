% This script changes the point cloud based on the uncertainty of the points, the angle (theta) from the scanner to the point and the laser radius at the hit location
% Only the range error along the scanner vector is simulated
% Diagnostics shows how the points are shifted, as well as the cylinder geometry if it is given

% Note: The cylinder centre, radius and direction are only used for the plot

function Shifted_Point_Cloud_Coord = Gaussian_Point_Sampling(Point_Cloud_Distributions, Plot)

    %% Inputs %%
        % The distributions
        distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
        distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
        distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;
        distribution_Sigma_cell     = Point_Cloud_Distributions.distribution_Sigma_cell;
        number_distributions_list   = Point_Cloud_Distributions.number_distributions_list;
        number_distributions        = Point_Cloud_Distributions.number_distributions;

    %% Sampling the distributions %%
        % Sampling the distributions
        Sampling_fun            = @(mu, Sigma) mvnrnd(mu, Sigma, 1);
        distribution_mu_cell_s  = cellfun(Sampling_fun, distribution_mu_cell, distribution_Sigma_cell, 'UniformOutput', false);
        mu_matrix_s             = vertcat(distribution_mu_cell_s{:});

        % Converted to the point cloud format
        number_scanners             = length(number_distributions_list);
        cumulative_num_distr_list   = [0, cumsum(number_distributions_list)];

        point_cloud_cell_s          = cell(1, number_scanners);

        for s = 1 : number_scanners
            first_ind   = cumulative_num_distr_list(s) + 1;
            end_ind     = cumulative_num_distr_list(s + 1);

            point_cloud_cell_s{s} = mu_matrix_s(first_ind : end_ind, :);
        end

        % A structure is created containing the shifted point cloud
        Shifted_Point_Cloud_Coord = struct('number_points_list', number_distributions_list, 'number_points_total', sum(number_distributions_list), 'point_cloud_cell', {point_cloud_cell_s});
        
    %% Plot %%
        if Plot == true
            % Scanner colour map
            scanner_cmap = cbrewer('qual', 'Set1', max(number_scanners, 3));
            scanner_cmap = max(scanner_cmap, 0);
            scanner_cmap = min(scanner_cmap, 1);
            
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
            set(gcf, 'color', [1, 1, 1])

            hold on
            grid on

            % The distributions
            for d = 1 : number_distributions
                % Distribution properties
                distribution_mu         = distribution_mu_cell{d};
                distribution_sigmae     = distribution_sigmae_cell{d};
                distribution_axes       = distribution_axes_cell{d};

                % Surface
                number_coord                            = 1e2;
                [ellipsoid_coord_matrix, number_coord]  = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_sigmae, distribution_axes, number_coord);

                x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.10, 'DisplayName', sprintf('Point cloud, 1 %s', '\sigma'));

                % Mu
                sc_mu = scatter3(distribution_mu(1), distribution_mu(2), distribution_mu(3), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k', 'DisplayName', '\mu');

                if d > 1
                    sc_mu.HandleVisibility      = 'Off';
                    surf_el.HandleVisibility    = 'Off';
                end
            end

            % The shifted point cloud
            for s = 1 : number_scanners
                scanner_colour = scanner_cmap(s, :);
                
                [x_s_points_list, y_s_points_list, z_s_points_list] = deal(point_cloud_cell_s{s}(:, 1), point_cloud_cell_s{s}(:, 2), point_cloud_cell_s{s}(:, 3));
                scatter3(x_s_points_list, y_s_points_list, z_s_points_list, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', scanner_colour, 'DisplayName', sprintf('Scanner %g, shifted', s));            
            end
            
            % Axes
            axis equal
            
            xlabel('x [m]');
            ylabel('y [m]');
            xlabel('z [m]');
            
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
            
            view(45, 45)

            % Legend
            legend('show', 'Location', 'Eastoutside');
            hold off

            disp('The point sampling script will finish (and close the plot) if you press a key on the keyboard')
            pause               

            close(1);
        end
end