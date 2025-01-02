% This script samples one point of the given multivariate normal distributions (3D ellipsoids)

function sampled_points_matrix = Random_Multivariate_Normal_Sampling(ellipsoid_centre_matrix, ellipsoid_sigmae_cell, ellipsoid_axes_cell, num_ellipsoids)

    %% Manual inputs %%
        Diagnostics = false;     % [true, false] Creates a diagnostic plot of the 

    %% Multivariate normal sampling %%
        % Conversion of centre matrix to cell
        num_dim                 = size(ellipsoid_centre_matrix, 2);
        ellipsoid_centre_cell   = mat2cell(ellipsoid_centre_matrix, ones(1, num_ellipsoids), num_dim);

        % Sampling over each ellipsoid
        number_samples  = 1;        
        Multivariate_Normal_Sampling_fun = @(ellipsoid_centre, ellipsoid_sigmae, ellipsoid_axes) Multivariate_Normal_Sampling(ellipsoid_centre, ellipsoid_sigmae, ellipsoid_axes);
        
        sampled_points_cell     = cellfun(Multivariate_Normal_Sampling_fun, ellipsoid_centre_cell, ellipsoid_sigmae_cell, ellipsoid_axes_cell, 'UniformOutput', false);
        sampled_points_matrix   = cell2mat(sampled_points_cell);
        
        % Local sampling function
        function sampled_point = Multivariate_Normal_Sampling(ellipsoid_centre, ellipsoid_sigmae, ellipsoid_axes)
            % Random distance along each axis
            covariance_matrix   = ellipsoid_sigmae.^2 .* eye(num_dim);
            origin              = zeros(1, num_dim);
            distance_vector     = mvnrnd(origin, covariance_matrix, number_samples);
            
            % Resulting sample
            sampled_point       = sum(distance_vector' .* ellipsoid_axes, 1) + ellipsoid_centre;
        end
        
    %% Diagnostics plot %%
        if Diagnostics == true
            number_coord    = 1e3;      % The number of points into which the ellipsoids are discretised
            num_STD         = 2;        % The number of standard deviations used for the ellipsoids in the plot
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
              
            % The point cloud and its uncertainty
            for e = 1: num_ellipsoids
                % Ellipsoid properties
                ellipsoid_centre    = ellipsoid_centre_matrix(e, :);
                ellipsoid_sigmae    = ellipsoid_sigmae_cell{e};
                ellipsoid_axes      = ellipsoid_axes_cell{e};
                
                % Surface
                [ellipsoid_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, num_STD*ellipsoid_sigmae, ellipsoid_axes, number_coord);

                x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', sprintf('Point cloud uncertainty, %g \\sigma', num_STD));

                if e > 1
                    surf_el.HandleVisibility    = 'Off';
                end
            end
            
            % Sampled points
            scatter3(sampled_points_matrix(:, 1), sampled_points_matrix(:, 2), sampled_points_matrix(:, 3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Sampled points');
            
            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')

            data_matrix = [ellipsoid_centre_matrix; sampled_points_matrix];
            data_max    = max(data_matrix, [], 1);
            data_min    = min(data_matrix, [], 1);
            data_ampl   = data_max - data_min;
            
            xlim([data_min(1), data_max(1)] + 0.1*data_ampl(1));
            ylim([data_min(2), data_max(2)] + 0.1*data_ampl(2));            
            zlim([data_min(3), data_max(3)] + 0.1*data_ampl(3));     
            
            % Aspect ratio
            axis equal

            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            disp('The script will continue and the plot will close when a key is pressed')
            pause()

            close(1)            
        end
        
end