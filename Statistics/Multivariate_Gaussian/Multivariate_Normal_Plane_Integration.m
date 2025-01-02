% This script determines the plane integrals of a 3D multivariate normal probability density function using geometric projection
% The plane lies on the point, and is normal to the given vector
% The normalised likelihood values are the density values, but made independent of dimension and uncertainty

function [f_plane_list, L_N_plane_list] = Multivariate_Normal_Plane_Integration(proj_points_matrix, alpha, proj_ellipse_centre, proj_ellipsoid_sigmae, proj_ellipsoid_axes, Plot)
        
    %% Coverage probability %%
        % The number of standard deviations follows from the coverage probability
        [num_points, num_dim]   = size(proj_points_matrix);
        P                       = 1 - alpha;
        m_STD                   = sqrt(chi2inv(P, num_dim));

    %% Integration onto the plane %%
        % The first two dimensions are used
        cross_section_matrix    = proj_points_matrix(:, 1 : num_dim - 1);
        mu_cross_section        = proj_ellipse_centre(1 : num_dim - 1);
        proj_ellipse_sigmae     = proj_ellipsoid_sigmae(1 : num_dim - 1);
        proj_ellipse_axes       = proj_ellipsoid_axes(:, 1 : num_dim - 1);
        
        % The covariance matrix is computed
        Sigma                   = proj_ellipse_sigmae.^2 .* eye(num_dim - 1);
        covariance_matrix       = proj_ellipse_axes \ Sigma * proj_ellipse_axes;

        % Normalised likelihood values
        delta_cell              = mat2cell(cross_section_matrix - mu_cross_section, ones(1, num_points), num_dim - 1);
        Mahalanobis_sq_fun      = @(delta) delta / covariance_matrix * delta';          % Note that this produces the Mahalanobis distance squared
        Mahalanobis_sq_list     = cellfun(Mahalanobis_sq_fun, delta_cell);
        
        L_N_plane_list          = exp(-1/2 * Mahalanobis_sq_list);
        
        % The probability density values
        normalisation_factor    = 1/(2*pi * prod(proj_ellipse_sigmae));
        dimension_factor        = P / chi2cdf(m_STD^2, 2);              % As the Gaussian is now considered 2-dimensional
        f_plane_list            = normalisation_factor * dimension_factor * L_N_plane_list;
        
    %% Plot %%
        if Plot == true
            % Message showing the integrated density value
            proj_ellipse_radii  = max(cross_section_matrix, [], 1);
            proj_ellipse_area   = pi * prod(proj_ellipse_radii);
            fprintf('The integrated density over the plane: %.3g \n', mean(f_plane_list) * proj_ellipse_area);
            
            % Colormap for the density values
            number_colours  = 1e3;
            density_cmap    = cbrewer('seq', 'Greens', number_colours);
            density_cmap    = max(density_cmap, 0);
            density_cmap    = min(density_cmap, 1);
            
            % The colour for each density value is found
            density_UB      = max(f_plane_list);
        
            colour_ind      = round(f_plane_list / density_UB * (number_colours - 1)) + 1;
            colour_matrix   = density_cmap(colour_ind, :);
            
            % Projected ellipse coordinates
            number_coord                = 1e3;
            proj_ellipse_coord_matrix   = Ellipse_Coordinate_Generator(proj_ellipse_centre, proj_ellipsoid_axes, m_STD*proj_ellipsoid_sigmae, number_coord); 

            % The bounds of the plot
            data_matrix     = [proj_ellipse_coord_matrix; proj_points_matrix];
            
            data_max_list   = max(data_matrix, [], 1); 
            data_min_list   = min(data_matrix, [], 1);

            x_lim = [data_min_list(1), data_max_list(1)] + 0.1*[-1, 1]*(data_max_list(1) - data_min_list(1));
            y_lim = [data_min_list(2), data_max_list(2)] + 0.1*[-1, 1]*(data_max_list(2) - data_min_list(2));

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            % First set of axes for geometry
            ax1 = axes;
            grid on
            hold on
            
            proj_el_pl_name = sprintf('Proj. ellipse, %.3g \\sigma', m_STD);
            proj_el_pl = plot(ax1, proj_ellipse_coord_matrix(:, 1), proj_ellipse_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'k', 'DisplayName', proj_el_pl_name);

            % The samples and their densities use a second set of axes
            ax2 = axes;
            
            samples_sc_name = 'Samples';
            samples_sc = scatter(ax2, cross_section_matrix(:, 1), cross_section_matrix(:, 2), 1e2, colour_matrix, 'filled', 'DisplayName', samples_sc_name);
            
            % The colourmap is applied to the second set of axes for the sample data
            colormap(ax2, density_cmap);
            cb = colorbar(ax2);
            caxis manual
            shading interp
            caxis([0, density_UB])
            colortitlehandle    = get(cb, 'Title');
            titlestring         = 'f [m^{-2}]';
            set(colortitlehandle, 'String', titlestring);
            cb.FontSize         = 15;
    
            % Text
            legend(ax1, [proj_el_pl, samples_sc], {proj_el_pl_name, samples_sc_name}, 'Location', 'Northoutside');

            xlabel(ax1, 'x [m]');
            ylabel(ax1, 'y [m]');

            set(ax1, 'FontSize', 15);
            set(ax1, 'LineWidth', 2);
        
            % Axes dimensions      
            xlim([ax1, ax2], x_lim);
            ylim([ax1, ax2], y_lim);

            ax1.DataAspectRatio = [1, 1, 1];
            ax2.DataAspectRatio = [1, 1, 1];

            % Change the size, link the axes and make the second set invisible
            axis_dimensions = [ax1.Position; ax2.Position];
            axis_starts     = max(axis_dimensions(:, 1:2), [], 1);
            axis_sizes      = min(axis_dimensions(:, 3:4), [], 1);    
            ax1.Position    = [axis_starts, axis_sizes];
            ax2.Position    = [axis_starts, axis_sizes];
            linkaxes([ax1, ax2]);
            ax2.Visible     = 'off';
        
            % Pause
            disp('The script will continue and the figures will close upon a button-press');
            pause();
            
            close(1);                
        end
        
end