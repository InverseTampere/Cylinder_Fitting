% This script computes the extrema of a 3D ellipsoid, w.r.t. the given vectors passing through the given point
% The extrema lie on the ellipsoid itself, the projections on the given vectors

function [extrema_matrix, proj_extrema_matrix, proj_ellipsoid_extents, proj_ellipsoid_centres] = Ellipsoid_Extrema(vector_matrix, point, ellipsoid_centre, ellipsoid_radii, ellipsoid_axes)

    %% Manual inputs %%
        Diagnostics = false;             % [true, false] Shows the ellipsoid and its extrema
        
    %% Extrema on the vectors %%
        % The extrema on the vector equal the projection of the ellipsoid on said vectors
        [proj_ellipsoid_centres, proj_extrema_matrix, proj_ellipsoid_extents] = Ellipsoid_to_Vector_Projection(vector_matrix, point, ellipsoid_centre, ellipsoid_radii, ellipsoid_axes);

    %% Extrema on the ellipsoid %%
        % Check that the dimensionality is correct
        [num_vectors, num_dim] = size(vector_matrix);
        
        if num_dim ~= 3
            error('The extrema can only be computed for a 3-dimensional ellipsoid');
        end

        % Each extrema is placed on the ellipsoid separately
        num_extrema     = size(proj_extrema_matrix, 3);
        extrema_matrix  = zeros(num_vectors, num_dim, num_extrema);
        
        for e = 1 : num_extrema
            % This extrema's matrix
            proj_extrema_matrix_half        = proj_extrema_matrix(:, :, e);
            
            % The intersections of the plane orthogonal to the vector with the ellipsoid
            [ellipse_extrema_half, ~, ~]    = Plane_Ellipsoid_Intersection(proj_extrema_matrix_half, vector_matrix, ellipsoid_centre, ellipsoid_radii, ellipsoid_axes);
            
            % Due to rounding errors in Cholesky decomposition, the extrema might not be found immediately
            delta = 1;
            
            while sum(isnan(ellipse_extrema_half)) > 0
                delta = delta + 1e-3;       % As such, the radii are slightly increased
                
                [ellipse_extrema_half, ~, ~] = Plane_Ellipsoid_Intersection(proj_extrema_matrix_half, vector_matrix, ellipsoid_centre, delta * ellipsoid_radii, ellipsoid_axes);
            end
            
            extrema_matrix(:, :, e) = ellipse_extrema_half;
        end

    %% Diagnostics plot %%
        if Diagnostics == true && num_dim == 3      % Only used if the given ellipsoid is 3D
            % Ellipsoid coordinates
            number_coord = 1e3;
            [ellipsoid_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, number_coord);

            % Used for plotting
            extrema_names           = {'Minimum', 'Maximum'};
            extrema_colours         = {'c', 'm'};
            extrema_vector_colours  = {'b', 'r'};
            
            vector_colours = cbrewer('qual', 'Set2', max(3, num_vectors));
            vector_colours = max(vector_colours, 0);
            vector_colours = min(vector_colours, 1);

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on

            % Ellipsoid
            x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
            y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
            z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

            surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'g', 'FaceAlpha', 0.1, 'DisplayName', 'Ellipsoid');

            % Vectors and their extrema
            for v = 1 : num_vectors
                % This vector
                vector_colour   = vector_colours(v, :);
                vector          = vector_matrix(v, :);
                proj_extent     = proj_ellipsoid_extents(v);
                proj_point      = proj_ellipsoid_centres(v, :);
                
                vector_string = sprintf('Vector %g', v);
                plot3([point(1), proj_point(1) + proj_extent/2*vector(1)*[-1, 1]], [point(2), proj_point(2) + proj_extent/2*vector(2)*[-1, 1]], [point(3), proj_point(3) + proj_extent/2*vector(3)*[-1, 1]], 'color', vector_colour, 'LineWidth', 2, 'DisplayName', vector_string);
                
                % Its extrema
                for e = 1 : num_extrema
                    extrema_colour      = extrema_colours{e};
                    proj_extrema_name   = sprintf('%s on vector', extrema_names{e});
                    sc_proj_extr        = scatter3(proj_extrema_matrix(v, 1, e), proj_extrema_matrix(v, 2, e), proj_extrema_matrix(v, 3, e), 'filled', 'MarkerFaceColor', extrema_colour, 'DisplayName', proj_extrema_name);

                    extrema_vector_colour = extrema_vector_colours{e};
                    extrema_name    = sprintf('%s on ellipsoid', extrema_names{e});
                    sc_extr         = scatter3(extrema_matrix(v, 1, e), extrema_matrix(v, 2, e), extrema_matrix(v, 3, e), 'filled', 'MarkerFaceColor', extrema_vector_colour, 'DisplayName', extrema_name);
                         
                    pl_extr = plot3([extrema_matrix(v, 1, e), proj_extrema_matrix(v, 1, e)], [extrema_matrix(v, 2, e), proj_extrema_matrix(v, 2, e)], [extrema_matrix(v, 3, e), proj_extrema_matrix(v, 3, e)], 'color', 'k', 'LineWidth', 1, 'LineStyle', '--');
                    pl_extr.HandleVisibility = 'Off';
                    
                    if v < num_vectors
                        sc_proj_extr.HandleVisibility   = 'Off';
                        sc_extr.HandleVisibility        = 'Off';
                    end
                end
            end

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            extrema_matrix_plot = [proj_extrema_matrix(:, :, 1); proj_extrema_matrix(:, :, 2)];
            data_matrix         = [ellipsoid_coord_matrix; point; extrema_matrix_plot];
            
            max_list    = max(data_matrix, [], 1);
            min_list    = min(data_matrix, [], 1);
            ampl_list   = max_list - min_list;

            xlim(0.2*ampl_list(1)*[-1, 1] + [min_list(1), max_list(1)]);
            ylim(0.2*ampl_list(2)*[-1, 1] + [min_list(2), max_list(2)]);
            zlim(0.2*ampl_list(3)*[-1, 1] + [min_list(3), max_list(3)]);

            axis equal      % Note that this changes the axis limits again

            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
           
            % Pause
            disp('The script will continue and the plots will close upon a button-press.');
            pause();
            
            close(1);
        end
        
end