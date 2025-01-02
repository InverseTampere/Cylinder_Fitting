% This script analytically determines the integrals of a 3D multivariate normal distribution onto a plane and vector over the given points

% The statistical output structure contains the density and normalised likelihood values on the plane and vector 
% The normalised likelihood values are the density values, but made independent of dimension and uncertainty

% Note that the distribution's expected value, axes and standard deviations are only used for diagnostics and may be empty otherwise

function Statistical_Outputs = Multivariate_Normal_Integration(points_matrix, plane_normal_vector, Projected_Distributions, distribution_Sigma, distribution_mu, distribution_axes, distribution_sigmae, Diagnostics)

    %% Transformation to the plane's coordinate frame %%
        % Normalisation of the normal vector
        plane_normal_vector = plane_normal_vector / norm(plane_normal_vector);
        
        % Rotation of the given points to the plane's coordinate frame
        num_dim     = length(plane_normal_vector);
        origin      = zeros(1, num_dim);
        
        [points_matrix_plane, plane_vector_basis, ~] = Vector_Based_Rotation(points_matrix, plane_normal_vector, origin);
        x_list  = points_matrix_plane(:, 1);
        y_list  = points_matrix_plane(:, 2);
        z_list  = points_matrix_plane(:, 3);
        
        % The precision matrix in the plane's coordinate frame
        Sigma   = plane_vector_basis * distribution_Sigma / plane_vector_basis;
        P       = inv(Sigma);

    %% Integration onto plane %%
        % Integration onto the plane is a single integral in the third dimension for the transformed distribution using the following set of coefficients    
        [alpha_z, beta_z, gamma_z, delta_z, epsilon_z, omega_z] = deal(Projected_Distributions.Plane.Coefficients.alpha, Projected_Distributions.Plane.Coefficients.beta, Projected_Distributions.Plane.Coefficients.gamma, Projected_Distributions.Plane.Coefficients.delta, Projected_Distributions.Plane.Coefficients.epsilon, Projected_Distributions.Plane.Coefficients.omega);
        
        % Normalised likelihood and density
        L_N_plane_list  = exp(alpha_z*x_list.^2 + beta_z*y_list.^2 + gamma_z*x_list.*y_list + delta_z*x_list + epsilon_z*y_list + omega_z);
        nu_plane        = 1/(2*pi * sqrt(P(3, 3) * det(Sigma)));      
        f_plane_list    = nu_plane * L_N_plane_list;
                
    %% Integration onto vector %%
        % Integration onto the vector is a double integral in the first and second dimension for the transformed distribution using the following set of coefficients    
        % The first integral is performed with respect to x and produces an ellipse on the y-z plane, the second integral uses coefficients which are a function of the first set forming a '1D ellipse' zeta*z^2 + eta*z + theta
        [alpha_x, zeta_y, eta_y, theta_y] = deal(Projected_Distributions.Vector.Coefficients.alpha, Projected_Distributions.Vector.Coefficients.zeta, Projected_Distributions.Vector.Coefficients.eta, Projected_Distributions.Vector.Coefficients.theta);
        
        % Normalised likelihood and density
        L_N_vector_list = exp(zeta_y*z_list.^2 + eta_y*z_list + theta_y);
        nu_vector       = 1/sqrt(-4*pi*alpha_x*P(1, 1)*det(Sigma));      
        f_vector_list   = nu_vector * L_N_vector_list;
        
    %% Output structures %%
        % Normalised likelihood and probability density values for the given points
        Statistical_Outputs = struct('Plane', struct('normalised_likelihood', L_N_plane_list, 'probability_density', f_plane_list), 'Vector', struct('normalised_likelihood', L_N_vector_list, 'probability_density', f_vector_list));
        
    %% Plot %%
        if Diagnostics == true
            %% Quantitative verification %%
                disp('The probability density for the points on the projected distributions should equal the values calculated above, i.e. MAE should be 0');
                
                %--% Plane %--%
                Proj_Diagnostics = false;
                [plane_proj_matrix, ~, height_list] = Point_to_Plane_Projection(points_matrix, plane_vector_basis, origin, Proj_Diagnostics);
    
                f_plane_list_3D     = sqrt(2*pi) * mvnpdf(plane_proj_matrix, mu_plane_3D, Sigma_plane_3D);      % Note that it has to be adjusted for dimensionality
                    
                MAE_plane           = mean(abs(f_plane_list - f_plane_list_3D));
                fprintf('MAE for the plane: %.3g \n', MAE_plane);

                %--% Vector %--%
                [~, mu_height]      = Point_to_Vector_Projection(mu_vector_3D, plane_normal_vector, origin);
                f_vector_list_3D    = normpdf(height_list, mu_height, sigma_vector);
                
                MAE_vector          = mean(abs(f_vector_list_3D - f_vector_list));
                fprintf('MAE for the vector: %.3g \n', MAE_vector);
                
            %% Visual inspection %%
                % Discretisation of the projected distributions and colour map
                n = 1e3;
                
                % First the vector plot is created, then the plane
                plot_names_cell = {'Vector', 'Plane'};
                num_plots       = length(plot_names_cell);
                
                % Number of standard deviations to cover all points
                m_STD_list      = [sqrt(chi2inv(1 - min(f_vector_list), 1)), sqrt(chi2inv(1 - min(f_plane_list), 2))];
                
                % Colour map for the density values
                density_cmap    = cbrewer('seq', 'Greens', n);
                density_cmap    = max(density_cmap, 0);
                density_cmap    = min(density_cmap, 1);
                
                % The colour bar limits
                f_maximum_list  = [max(f_vector_list), max(f_plane_list)];
                
                %--% Original coordinate frame %--%
                for p = 1 : num_plots
                    % This plot's information
                    plot_name   = plot_names_cell{p};
                    m_STD       = m_STD_list(p);
                    f_maximum   = f_maximum_list(p);
                    
                    figure(p)
                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                    set(gcf, 'color', [1, 1, 1])    
                    
                    sgtitle(plot_name);

                    % First set of axes for geometry
                    ax1 = axes;
                    hold on
                    grid on

                    % The distribution
                    distribution_radii = distribution_sigmae * m_STD;
                    [distr_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, n);

                    a_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                    b_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                    c_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                    surf_distr_name = sprintf('Distribution, %.2g %s', m_STD, '\sigma');
                    surf_distr      = surf(ax1, a_distribution, b_distribution, c_distribution, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', surf_distr_name);

                    % The plane
                    plane_corner_matrix     = Plane_Corner_Points(plane_normal_vector, origin, [distr_coord_matrix; points_matrix]);
                    pa_plane_name           = 'Plane';
                    pa_plane                = patch(ax1, plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'g', 'FaceAlpha', 0.25, 'DisplayName', pa_plane_name);   

                    % The vector
                    vector_length   = max(abs([height_list, mu_height]));
                    pl_vector_name  = 'Normal vector';
                    pl_vector       = plot3(vector_length * [-1, 1] * plane_normal_vector(1), vector_length * [-1, 1] * plane_normal_vector(2), vector_length * [-1, 1] * plane_normal_vector(3), 'LineWidth', 2, 'color', 'm', 'DisplayName', pl_vector_name);

                    % The points
                    sc_points_name  = 'Points';
                    sc_points       = scatter3(ax1, points_matrix(:, 1), points_matrix(:, 2), points_matrix(:, 3), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', sc_points_name);

                    % Second set of axes for the density values
                    ax2 = axes;
                    hold on
                    
                    if strcmp(plot_name, 'Plane')
                        % Points projected onto the plane
                        density_ind         = round(f_plane_list / f_maximum * (n - 1)) + 1;
                        density_colours     = density_cmap(density_ind, :);

                        sc_pl_points_name   = 'Plane projection';
                        sc_pl_points        = scatter3(plane_proj_matrix(:, 1), plane_proj_matrix(:, 2), plane_proj_matrix(:, 3), 1e2, density_colours, 'filled', 'MarkerEdgeColor', 'r', 'DisplayName', sc_pl_points_name);

                        % The colour map is applied to the second set of axes
                        colormap(ax2, density_cmap);
                        cb = colorbar(ax2);
                        caxis manual
                        shading interp
                        caxis([0, f_maximum]);
                        colortitlehandle    = get(cb, 'Title');
                        titlestring         = 'f [m^{-2}]';
                        set(colortitlehandle, 'String', titlestring);
                        cb.FontSize         = 15;

                        % Text
                        legend(ax1, [surf_distr, pa_plane, pl_vector, sc_points, sc_pl_points], {surf_distr_name, pa_plane_name, pl_vector_name, sc_points_name, sc_pl_points_name}, 'Location', 'Northoutside');
                    
                    elseif strcmp(plot_name, 'Vector')
                        % Projected onto the vector
                        density_ind         = round(f_vector_list / f_maximum * (n - 1)) + 1;
                        density_colours     = density_cmap(density_ind, :);

                        vector_proj_matrix  = height_list .* plane_normal_vector;
                        sc_ve_points_name   = 'Vector projection';
                        sc_ve_points        = scatter3(vector_proj_matrix(:, 1), vector_proj_matrix(:, 2), vector_proj_matrix(:, 3), 1e2, density_colours, 'filled', 'MarkerEdgeColor', 'c', 'DisplayName', sc_ve_points_name);
                        
                        % The colour map is applied to the second set of axes
                        colormap(ax2, density_cmap);
                        cb = colorbar(ax2);
                        caxis manual
                        shading interp
                        caxis([0, f_maximum]);
                        colortitlehandle    = get(cb, 'Title');
                        titlestring         = 'f [1/m]';
                        set(colortitlehandle, 'String', titlestring);
                        cb.FontSize         = 15;

                        % Text
                        legend(ax1, [surf_distr, pa_plane, pl_vector, sc_points, sc_ve_points], {surf_distr_name, pa_plane_name, pl_vector_name, sc_points_name, sc_ve_points_name}, 'Location', 'Northoutside');
                    end
                    
                    xlabel(ax1, 'a [m]');
                    ylabel(ax1, 'b [m]');
                    zlabel(ax1, 'c [m]');

                    set(ax1, 'FontSize', 15);
                    set(ax1, 'LineWidth', 2);

                    % Axes dimensions
                    ax1.DataAspectRatio = [1, 1, 1];
                    ax2.DataAspectRatio = [1, 1, 1];

                    x_lim                   = [min([ax1.XLim, ax2.XLim]), max([ax1.XLim, ax2.XLim])];
                    y_lim                   = [min([ax1.YLim, ax2.YLim]), max([ax1.YLim, ax2.YLim])];
                    z_lim                   = [min([ax1.ZLim, ax2.ZLim]), max([ax1.ZLim, ax2.ZLim])];

                    [ax1.XLim, ax2.XLim]    = deal(x_lim);
                    [ax1.YLim, ax2.YLim]    = deal(y_lim);
                    [ax1.ZLim, ax2.ZLim]    = deal(z_lim);

                    % Change the size, link the axes and make the second set invisible                
                    axis_dimensions = [ax1.Position; ax2.Position];
                    axis_starts     = max(axis_dimensions(:, 1:2), [], 1);
                    axis_sizes      = min(axis_dimensions(:, 3:4), [], 1);    
                    ax1.Position    = [axis_starts, axis_sizes];
                    ax2.Position    = [axis_starts, axis_sizes];
                    linkaxes([ax1, ax2]);
                    ax2.Visible     = 'off';

                    % The perspective
                    linkprop([ax1, ax2], 'view');       % Links the perspective of the two axes
                    view(ax1, 45, 45);
                end
                
                %--% Plane-aligned coordinate frame %--%
                for p = 1 : num_plots
                    % This plot's information
                    plot_name   = plot_names_cell{p};
                    m_STD       = m_STD_list(p);
                    f_maximum   = f_maximum_list(p);
                    
                    figure(p + 2)
                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                    set(gcf, 'color', [1, 1, 1])    
                    
                    sgtitle(plot_name);

                    % First set of axes for geometry
                    ax1 = axes;
                    hold on
                    grid on

                    % The distribution
                    mu                      = (plane_vector_basis * distribution_mu')';
        
                    distribution_radii      = distribution_sigmae * m_STD;
                    distribution_axes_xyz   = (plane_vector_basis * distribution_axes')';
                    [distr_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(mu, distribution_radii, distribution_axes_xyz, n);

                    x_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                    y_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                    z_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                    surf_distr_name = sprintf('Distribution, %.2g %s', m_STD, '\sigma');
                    surf_distr      = surf(ax1, x_distribution, y_distribution, z_distribution, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', surf_distr_name);

                    % The points
                    sc_points_name  = 'Points';
                    sc_points       = scatter3(ax1, x_list, y_list, z_list, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', sc_points_name);

                    % Projected plane distributions
                    mu_plane            = Projected_Distributions.Plane.Projection.mu;
                    distr_axes_plane    = Projected_Distributions.Plane.Projection.distr_axes;
                    sigmae_plane        = Projected_Distributions.Plane.Projection.sigmae;
                    radii_plane         = m_STD * sigmae_plane;
                    
                    plane_distr_coord   = Ellipse_Coordinate_Generator(mu_plane, distr_axes_plane, radii_plane, n);

                    pl_plane_name   = sprintf('Plane distribution, %.2g %s', m_STD, '\sigma');
                    pl_plane        = plot3(plane_distr_coord(:, 1), plane_distr_coord(:, 2), zeros(n, 1), 'LineWidth', 2, 'color', 'g', 'DisplayName', pl_plane_name);

                    % Vector distribution expected value
                    mu_vector       = Projected_Distributions.Vector.Projection.mu;
                    
                    sc_mu_vec_name  = '\mu_{v}';
                    sc_mu_vec       = scatter3(ax1, 0, 0, mu_vector, 'filled', 'MarkerFaceColor', 'm', 'DisplayName', sc_mu_vec_name);

                    % Second set of axes for the density values
                    ax2 = axes;
                    hold on
                    
                    if strcmp(plot_name, 'Plane')
                        % Points projected onto the plane
                        density_ind         = round(f_plane_list / f_maximum * (n - 1)) + 1;
                        density_colours     = density_cmap(density_ind, :);

                        sc_pl_points_name   = 'Plane projection';
                        num_points          = size(points_matrix, 1);
                        sc_pl_points        = scatter3(x_list, y_list, zeros(num_points, 1), 1e2, density_colours, 'filled', 'MarkerEdgeColor', 'r', 'DisplayName', sc_pl_points_name);

                        % The colour map is applied to the second set of axes
                        colormap(ax2, density_cmap);
                        cb = colorbar(ax2);
                        caxis manual
                        shading interp
                        caxis([0, f_maximum]);
                        colortitlehandle    = get(cb, 'Title');
                        titlestring         = 'f [m^{-2}]';
                        set(colortitlehandle, 'String', titlestring);
                        cb.FontSize         = 15;

                        % Text
                        legend(ax1, [surf_distr, pa_plane, pl_vector, sc_points, pl_plane, sc_mu_vec, sc_pl_points], {surf_distr_name, pa_plane_name, pl_vector_name, sc_points_name, pl_plane_name, sc_mu_vec_name, sc_pl_points_name}, 'Location', 'Northoutside');
                    
                    elseif strcmp(plot_name, 'Vector')
                        % Projected onto the vector
                        density_ind         = round(f_vector_list / f_maximum * (n - 1)) + 1;
                        density_colours     = density_cmap(density_ind, :);

                        sc_ve_points_name   = 'Vector projection';
                        num_points          = size(points_matrix, 1);
                        sc_ve_points        = scatter3(zeros(num_points, 1), zeros(num_points, 1), z_list, 1e2, density_colours, 'filled', 'MarkerEdgeColor', 'c', 'DisplayName', sc_ve_points_name);
                        
                        % The colour map is applied to the second set of axes
                        colormap(ax2, density_cmap);
                        cb = colorbar(ax2);
                        caxis manual
                        shading interp
                        caxis([0, f_maximum]);
                        colortitlehandle    = get(cb, 'Title');
                        titlestring         = 'f [1/m]';
                        set(colortitlehandle, 'String', titlestring);
                        cb.FontSize         = 15;

                        % Text
                        legend(ax1, [surf_distr, sc_points, pl_plane, sc_mu_vec, sc_ve_points], {surf_distr_name, sc_points_name, pl_plane_name, sc_mu_vec_name, sc_ve_points_name}, 'Location', 'Northoutside');
                    end
                    
                    xlabel(ax1, 'x [m]');
                    ylabel(ax1, 'y [m]');
                    zlabel(ax1, 'z [m]');

                    set(ax1, 'FontSize', 15);
                    set(ax1, 'LineWidth', 2);

                    % Axes dimensions
                    ax1.DataAspectRatio = [1, 1, 1];
                    ax2.DataAspectRatio = [1, 1, 1];

                    x_lim                   = [min([ax1.XLim, ax2.XLim]), max([ax1.XLim, ax2.XLim])];
                    y_lim                   = [min([ax1.YLim, ax2.YLim]), max([ax1.YLim, ax2.YLim])];
                    z_lim                   = [min([ax1.ZLim, ax2.ZLim]), max([ax1.ZLim, ax2.ZLim])];

                    [ax1.XLim, ax2.XLim]    = deal(x_lim);
                    [ax1.YLim, ax2.YLim]    = deal(y_lim);
                    [ax1.ZLim, ax2.ZLim]    = deal(z_lim);

                    % Change the size, link the axes and make the second set invisible                
                    axis_dimensions = [ax1.Position; ax2.Position];
                    axis_starts     = max(axis_dimensions(:, 1:2), [], 1);
                    axis_sizes      = min(axis_dimensions(:, 3:4), [], 1);    
                    ax1.Position    = [axis_starts, axis_sizes];
                    ax2.Position    = [axis_starts, axis_sizes];
                    linkaxes([ax1, ax2]);
                    ax2.Visible     = 'off';

                    % The perspective
                    linkprop([ax1, ax2], 'view');       % Links the perspective of the two axes
                    view(ax1, 45, 45);
                end

                %--% Pause %--%
                disp('The script will continue and the figures will close upon a button-press');
                pause();

                close(1:4);                
        end
        
end