% The projection of the ellipsoids onto a plane is the maximum area of its intersection with said plane, where it can move along its normal vector

% Note that the coordinates are returned in the original coordinate frame, and defined on the plane itself with the plane point as origin

function Projected_Ellipsoids = Ellipsoid_to_Plane_Projection(ellipsoid_centre_cell, ellipsoid_axes_cell, ellipsoid_radii_cell, plane_point, plane_vector_basis, Plot)

    %% Manual inputs %%        
        Diagnostics  = false;        % [true, false] Shows the projected points onto the plane, used to fit the projected ellipse, and the end result for an individual ellipse
        
    %% Ellipsoid projection %%
        % Results are provided both in the local (plane point as origin, plane normal vector as z-axis) and global frames
        Ellipsoid_Projection_fun = @(ellipsoid_centre, ellipsoid_axes, ellipsoid_radii) Ellipsoid_Projection(ellipsoid_centre, ellipsoid_axes, ellipsoid_radii, plane_point, plane_vector_basis, Diagnostics);

        [proj_ellipse_radii_cell, proj_ellipse_centre_cell, proj_ellipse_axes_cell, proj_ellipse_centre_3D_cell, proj_ellipse_axes_3D_cell] = cellfun(Ellipsoid_Projection_fun, ellipsoid_centre_cell, ellipsoid_axes_cell, ellipsoid_radii_cell, 'UniformOutput', false);
       
        % Results are stored in a structure
        number_ellipsoids       = length(ellipsoid_centre_cell);
        Projected_Ellipsoids    = struct('number_ellipsoids', number_ellipsoids, 'plane_point', plane_point, 'projection_vector_basis', plane_vector_basis, ...
                                         'proj_ellipse_radii_cell', {proj_ellipse_radii_cell}, 'proj_ellipse_centre_cell', {proj_ellipse_centre_cell}, 'proj_ellipse_axes_cell', {proj_ellipse_axes_cell}, 'proj_ellipse_centre_3D_cell', {proj_ellipse_centre_3D_cell}, 'proj_ellipse_axes_3D_cell', {proj_ellipse_axes_3D_cell});

    %% Plot %%
        % All ellipsoids and their projections are shown
        if Plot == true
            num_coord   = 1e3;
            
            % The centres and radii of the (projected) ellipsoids
            ellipsoid_centre_matrix         = vertcat(ellipsoid_centre_cell{:});
            proj_ellipse_centre_matrix      = vertcat(proj_ellipse_centre_cell{:});
            proj_ellipse_radii_matrix       = vertcat(proj_ellipse_radii_cell{:});
            
            ellipsoid_centre_matrix_r       = (plane_vector_basis * (ellipsoid_centre_matrix - plane_point)')';
            proj_ellipse_centre_3D_matrix   = vertcat(proj_ellipse_centre_3D_cell{:});
            
            % The data bounds
            max_radius  = max(proj_ellipse_radii_matrix, [], 'all');
            data_LB     = min([ellipsoid_centre_matrix; proj_ellipse_centre_3D_matrix], [], 1) - max_radius;
            data_UB     = max([ellipsoid_centre_matrix; proj_ellipse_centre_3D_matrix], [], 1) + max_radius;
            
            proj_data_LB = min([ellipsoid_centre_matrix_r(:, 1:num_dim - 1); proj_ellipse_centre_matrix], [], 1) - max_radius;
            proj_data_UB = max([ellipsoid_centre_matrix_r(:, 1:num_dim - 1); proj_ellipse_centre_matrix], [], 1) + max_radius;
            
            % Plane coordinates
            num_dim                 = length(plane_point);
            plane_normal_vector     = plane_vector_basis(num_dim, :);
            plane_corner_matrix_3D  = Plane_Corner_Points(plane_normal_vector, plane_point, [data_LB; data_UB]);

            origin                      = zeros(1, num_dim);
            z_axis                      = [0, 0, 1];
            proj_plane_corner_matrix    = Plane_Corner_Points(z_axis, origin, [proj_data_LB, 0; proj_data_UB, 0]);  

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   
            
            %--% Geometry %--%
            for e = 1 : number_ellipsoids
                % Ellipsoid coordinates    
                ellipsoid_centre    = ellipsoid_centre_cell{e};
                ellipsoid_radii     = ellipsoid_radii_cell{e};
                ellipsoid_axes      = ellipsoid_axes_cell{e};
                
                [ellipsoid_coord_matrix, num_coord] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, num_coord);
                ellipsoid_coord_matrix_t            = ellipsoid_coord_matrix - plane_point;
                ellipsoid_coord_matrix_r            = (plane_vector_basis * ellipsoid_coord_matrix_t')';    % In the plane's coordinate basis

                % Projected ellipse coordinates
                proj_ellipse_centre     = proj_ellipse_centre_cell{e};
                proj_ellipse_centre_3D  = proj_ellipse_centre_3D_cell{e};
                proj_ellipse_axes       = proj_ellipse_axes_cell{e};
                proj_ellipse_axes_3D    = proj_ellipse_axes_3D_cell{e};
                proj_ellipse_radii      = proj_ellipse_radii_cell{e};
                
                proj_ellipse_coordinates    = Ellipse_Coordinate_Generator(proj_ellipse_centre, proj_ellipse_axes, proj_ellipse_radii, num_coord);
                proj_ellipse_coordinates_3D = Ellipse_Coordinate_Generator(proj_ellipse_centre_3D, proj_ellipse_axes_3D, proj_ellipse_radii, num_coord);

                %--% Original coordinate frame %--%
                subplot(1, 2, 1)
                hold on
                grid on

                % Ellipsoid
                x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(num_coord) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(num_coord) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(num_coord) * [1, 1]);

                su_ell = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', 'Ellipsoid');

                % Projected ellipse
                pl_proj_ell = plot3(proj_ellipse_coordinates_3D(:, 1), proj_ellipse_coordinates_3D(:, 2), proj_ellipse_coordinates_3D(:, 3), 'LineWidth', 2, 'color', 'c', 'DisplayName', 'Projected ellipsoid');

                % Plane
                if e == 1
                    % Normal vector and point
                    scatter3(plane_point(1), plane_point(2), plane_point(3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Plane point');
                    
                    normal_vector_scaled = plane_normal_vector * sqrt(sum((data_UB - data_LB).^2));
                    pl_normal = plot3(plane_point(1) + normal_vector_scaled(1)*[-1, 1], plane_point(2) + normal_vector_scaled(2)*[-1, 1], plane_point(3) + normal_vector_scaled(3)*[-1, 1], 'color', 'r', 'LineWidth', 2, 'DisplayName', 'Normal vector');

                    % Plane
                    patch(plane_corner_matrix_3D(:, 1), plane_corner_matrix_3D(:, 2), plane_corner_matrix_3D(:, 3), 'r', 'FaceAlpha', 0.25, 'DisplayName', 'Plane');   
                end
                
                if e > 1
                    su_ell.HandleVisibility         = 'Off';
                    pl_proj_ell.HandleVisibility    = 'Off';
                    pl_normal.HandleVisibility      = 'Off';
                end
                
                %--% Plane's coordinate basis %--%
                subplot(1, 2, 2)
                hold on
                grid on

                % Ellipsoid
                x_ellipsoid = reshape(ellipsoid_coord_matrix_r(:, 1), sqrt(num_coord) * [1, 1]);
                y_ellipsoid = reshape(ellipsoid_coord_matrix_r(:, 2), sqrt(num_coord) * [1, 1]);
                z_ellipsoid = reshape(ellipsoid_coord_matrix_r(:, 3), sqrt(num_coord) * [1, 1]);

                surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', 'Ellipsoid');

                % The ellipsoid axes
                ellipsoid_axes_r    = (plane_vector_basis * ellipsoid_axes')';
                ellipsoid_centre_t  = ellipsoid_centre - plane_point;
                ellipsoid_centre_r  = (plane_vector_basis * ellipsoid_centre_t')';       

                % Projected ellipse
                plot3(proj_ellipse_coordinates(:, 1), proj_ellipse_coordinates(:, 2), zeros(num_coord, 1), 'LineWidth', 2, 'color', 'c', 'DisplayName', 'Projected ellipsoid');

                % Plane
                if e == 1
                    % Normal vector and point
                    scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Plane point');
                    
                    normal_vector_scaled = z_axis * sqrt(sum((proj_data_UB - proj_data_LB).^2));
                    pl_normal = plot3(normal_vector_scaled(1)*[-1, 1], normal_vector_scaled(2)*[-1, 1], normal_vector_scaled(3)*[-1, 1], 'color', 'r', 'LineWidth', 2, 'DisplayName', 'Normal vector');

                    % Plane
                    patch(proj_plane_corner_matrix(:, 1), proj_plane_corner_matrix(:, 2), proj_plane_corner_matrix(:, 3), 'r', 'FaceAlpha', 0.25, 'DisplayName', 'Plane');   
                end
            end
            
            %--% Plot formatting %--%
            subplot(1, 2, 1)
            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            ampl_list = data_UB - data_LB;
            
            xlim(0.2*ampl_list(1)*[-1, 1] + [data_LB(1), data_UB(1)]);
            ylim(0.2*ampl_list(2)*[-1, 1] + [data_LB(2), data_UB(2)]);
            zlim(0.2*ampl_list(3)*[-1, 1] + [data_LB(3), data_UB(3)]);

            axis equal      

            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
                
            subplot(1, 2, 2)
            % Axes
            xlabel('a [m]');
            ylabel('b [m]');
            zlabel('c [m]');

            axis equal    
            
            ampl_list = proj_data_UB - proj_data_LB;

            xlim(0.2*ampl_list(1)*[-1, 1] + [proj_data_LB(1), proj_data_UB(1)]);
            ylim(0.2*ampl_list(2)*[-1, 1] + [proj_data_LB(2), proj_data_UB(2)]);

            view(0, 90);

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Pause
            disp('The script will continue and the figure will close upon a button-press');
            pause();

            close(1);   
        end
                                  
    %% Local function that projects each ellipsoid onto the given plane %%
        % Outputs are given in the projected coordinate frame (2D) and in the original one (3D)
        function [proj_ellipse_radii, proj_ellipse_centre, proj_ellipse_axes, proj_ellipse_centre_3D, proj_ellipse_axes_3D] = Ellipsoid_Projection(ellipsoid_centre, ellipsoid_axes, ellipsoid_radii, plane_point, plane_vector_basis, Diagnostics)                    
                %--% Extrema %--%
                % Ellipsoid axes are normalised
                ellipsoid_axes = ellipsoid_axes ./ sqrt(sum(ellipsoid_axes.^2, 2));

                % The extrema on the ellipsoid w.r.t. the given vectors, in original coordinate frame
                [extrema_loc_matrix_3D, ~, ~]   = Ellipsoid_Extrema(plane_vector_basis, ellipsoid_centre, ellipsoid_centre, ellipsoid_radii, ellipsoid_axes);            
                extrema_loc_matrix_min          = extrema_loc_matrix_3D(:, :, 1);
                extrema_loc_matrix_max          = extrema_loc_matrix_3D(:, :, 2);
                
                % Transformed to the plane's coordinate basis
                extrema_loc_matrix_min_r = (plane_vector_basis * (extrema_loc_matrix_min - plane_point)')';
                extrema_loc_matrix_max_r = (plane_vector_basis * (extrema_loc_matrix_max - plane_point)')';
                
                extrema_loc_matrix = cat(3, extrema_loc_matrix_min_r, extrema_loc_matrix_max_r);
                    
                %--% Least-squares ellipse fitting %--%
                % For least-squares fitting one more point is needed on the plane, which is taken in-between the other two
                psi                         = (1 + sqrt(5))/2;      % The golden-ratio is used, rather arbitrarily, but it prevents singularities when the vectors are basis vectors
                additional_vector           = 1/psi * plane_vector_basis(1, :) + (1 - 1/psi) * plane_vector_basis(2, :);
                additional_points_matrix    = Ellipsoid_Extrema(additional_vector, ellipsoid_centre, ellipsoid_centre, ellipsoid_radii, ellipsoid_axes);
                
                % The full matrix of points, projected onto the plane
                num_dim                             = length(plane_point);
                ellipse_points_matrix               = [extrema_loc_matrix_min(1:num_dim - 1, :); extrema_loc_matrix_max(1:num_dim - 1, :); additional_points_matrix(:, :, 1)];            
                [~, ellipse_points_matrix_2D, ~]    = Point_to_Plane_Projection(ellipse_points_matrix, plane_vector_basis, plane_point, Diagnostics);
                
                % Least-squares fitting, which is exact in this case
                [proj_ellipse_centre, proj_ellipse_radii, proj_ellipse_axes] = Ellipse_Least_Squares_Fitting(ellipse_points_matrix_2D);
    
                % The centre and axes are transformed to the original coordinate frame
                proj_ellipse_centre_t   = [proj_ellipse_centre, 0];               % The third dimension is added again
                proj_ellipse_axes_r     = [proj_ellipse_axes, zeros(2, 1)];
                
                proj_ellipse_centre_3D  = (plane_vector_basis \ proj_ellipse_centre_t')' + plane_point;
                proj_ellipse_axes_3D    = (plane_vector_basis \ proj_ellipse_axes_r')';
                            
                %--% Diagnostics %--%
                if Diagnostics == true
                    % Ellipsoid coordinates    
                    num_coord = 1e3;
                    [ellipsoid_coord_matrix, num_coord] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, num_coord);         % Original coordinate frame
                    ellipsoid_coord_matrix_t            = ellipsoid_coord_matrix - plane_point;
                    ellipsoid_coord_matrix_r            = (plane_vector_basis * ellipsoid_coord_matrix_t')';                                                    % Projection coordinate frame
    
                    % Projected ellipse coordinates
                    proj_ellipse_coordinates_3D = Ellipse_Coordinate_Generator(proj_ellipse_centre_3D, proj_ellipse_axes_3D, proj_ellipse_radii, num_coord);
                    proj_ellipse_coordinates    = Ellipse_Coordinate_Generator(proj_ellipse_centre, proj_ellipse_axes, proj_ellipse_radii, num_coord);
    
                    % Plane
                    plane_normal_vector         = plane_vector_basis(num_dim, :);
                    plane_corner_matrix_3D      = Plane_Corner_Points(plane_normal_vector, plane_point, [ellipsoid_coord_matrix; proj_ellipse_coordinates_3D]);
                    origin                      = zeros(1, num_dim);
                    z_axis                      = [0, 0, 1];
                    proj_plane_corner_matrix    = Plane_Corner_Points(z_axis, origin, [ellipsoid_coord_matrix_r; proj_ellipse_coordinates, zeros(num_coord, 1)]);
    
                    figure(1)
                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                    set(gcf, 'color', [1, 1, 1])     
    
                    %--% Original coordinate system %--%
                    subplot(1, 2, 1)
                    hold on
                    grid on
    
                    % Ellipsoid
                    x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(num_coord) * [1, 1]);
                    y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(num_coord) * [1, 1]);
                    z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(num_coord) * [1, 1]);
    
                    surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', 'Ellipsoid');
    
                    % The ellipsoid axes
                    for d = 1 : num_dim
                        ellipsoid_axis = ellipsoid_radii(d) * ellipsoid_axes(d, :);   
                        e_ax = plot3(ellipsoid_centre(1) + [0, ellipsoid_axis(1)], ellipsoid_centre(2) + [0, ellipsoid_axis(2)], ellipsoid_centre(3) + [0, ellipsoid_axis(3)], 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Ellipsoid axes');
    
                        if d > 1
                            e_ax.HandleVisibility = 'Off';
                        end
                    end            
    
                    % Extrema
                    extrema_loc_matrix_plot = [extrema_loc_matrix_3D(:, :, 1); extrema_loc_matrix_3D(:, :, 2)];
                    scatter3(extrema_loc_matrix_plot(:, 1), extrema_loc_matrix_plot(:, 2), extrema_loc_matrix_plot(:, 3), 'filled', 'MarkerFaceColor', 'c', 'DisplayName', 'Extrema');
    
                    % Samples
                    scatter3(ellipse_points_matrix(:, 1), ellipse_points_matrix(:, 2), ellipse_points_matrix(:, 3), 'filled', 'MarkerFaceColor', 'm', 'DisplayName', 'Samples');
    
                    % Projected ellipse
                    plot3(proj_ellipse_coordinates_3D(:, 1), proj_ellipse_coordinates_3D(:, 2), proj_ellipse_coordinates_3D(:, 3), 'LineWidth', 2, 'color', 'c', 'DisplayName', 'Projected ellipsoid');
    
                    % Normal vector
                    normal_vector_scaled = plane_normal_vector * max(ellipsoid_radii);
                    plot3(ellipsoid_centre(1) + normal_vector_scaled(1)*[-1, 1], ellipsoid_centre(2) + normal_vector_scaled(2)*[-1, 1], ellipsoid_centre(3) + normal_vector_scaled(3)*[-1, 1], 'color', 'r', 'LineWidth', 2, 'DisplayName', 'Normal vector');
    
                    % Plane
                    scatter3(plane_point(1), plane_point(2), plane_point(3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Plane point');
                    patch(plane_corner_matrix_3D(:, 1), plane_corner_matrix_3D(:, 2), plane_corner_matrix_3D(:, 3), 'r', 'FaceAlpha', 0.25, 'DisplayName', 'Plane');   
    
                    % Axes
                    xlabel('x [m]');
                    ylabel('y [m]');
                    zlabel('z [m]');
    
                    max_list    = max([ellipsoid_coord_matrix; proj_ellipse_coordinates_3D], [], 1);
                    min_list    = min([ellipsoid_coord_matrix; proj_ellipse_coordinates_3D], [], 1);
                    ampl_list   = max_list - min_list;
    
                    xlim(0.2*ampl_list(1)*[-1, 1] + [min_list(1), max_list(1)]);
                    ylim(0.2*ampl_list(2)*[-1, 1] + [min_list(2), max_list(2)]);
                    zlim(0.2*ampl_list(3)*[-1, 1] + [min_list(3), max_list(3)]);
    
                    axis equal     

                    view(45, 45);
    
                    % Legend
                    legend('show', 'location', 'northoutside');
    
                    set(gca, 'FontSize', 15);
                    set(gca, 'LineWidth', 2);
    
                    hold off    
    
                    %--% Projection onto plane %--%
                    subplot(1, 2, 2)
                    hold on
                    grid on
    
                    % Ellipsoid
                    x_ellipsoid = reshape(ellipsoid_coord_matrix_r(:, 1), sqrt(num_coord) * [1, 1]);
                    y_ellipsoid = reshape(ellipsoid_coord_matrix_r(:, 2), sqrt(num_coord) * [1, 1]);
                    z_ellipsoid = reshape(ellipsoid_coord_matrix_r(:, 3), sqrt(num_coord) * [1, 1]);
    
                    surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'DisplayName', 'Ellipsoid');
    
                    % The ellipsoid axes
                    ellipsoid_axes_r    = (plane_vector_basis * ellipsoid_axes')';
                    ellipsoid_centre_t  = ellipsoid_centre - plane_point;
                    ellipsoid_centre_r  = (plane_vector_basis * ellipsoid_centre_t')';
                    
                    for d = 1 : num_dim
                        ellipsoid_axis = ellipsoid_radii(d) * ellipsoid_axes_r(d, :);   
                        e_ax = plot3(ellipsoid_centre_r(1) + [0, ellipsoid_axis(1)], ellipsoid_centre_r(2) + [0, ellipsoid_axis(2)], ellipsoid_centre_r(3) + [0, ellipsoid_axis(3)], 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Ellipsoid axes');
    
                        if d > 1
                            e_ax.HandleVisibility = 'Off';
                        end
                    end            
    
                    % Extrema
                    extrema_loc_matrix_plot = [extrema_loc_matrix(:, :, 1), extrema_loc_matrix(:, :, 2)];
                    scatter3(extrema_loc_matrix_plot(:, 1), extrema_loc_matrix_plot(:, 2), extrema_loc_matrix_plot(:, 3), 'filled', 'MarkerFaceColor', 'c', 'DisplayName', 'Extrema');
    
                    % Samples
                    scatter3(ellipse_points_matrix_2D(:, 1), ellipse_points_matrix_2D(:, 2), ellipse_points_matrix_2D(:, 3), 'filled', 'MarkerFaceColor', 'm', 'DisplayName', 'Samples');
    
                    % Projected ellipse
                    plot3(proj_ellipse_coordinates(:, 1), proj_ellipse_coordinates(:, 2), zeros(num_coord, 1), 'LineWidth', 2, 'color', 'c', 'DisplayName', 'Projected ellipsoid');
    
                    % Normal vector
                    normal_vector_scaled = z_axis * max(ellipsoid_radii);
                    plot3(ellipsoid_centre_r(1) + normal_vector_scaled(1)*[-1, 1], ellipsoid_centre_r(2) + normal_vector_scaled(2)*[-1, 1], ellipsoid_centre_r(3) + normal_vector_scaled(3)*[-1, 1], 'color', 'r', 'LineWidth', 2, 'DisplayName', 'Normal vector');
    
                    % Plane
                    scatter3(0, 0, 0, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Plane point');
                    patch(proj_plane_corner_matrix(:, 1), proj_plane_corner_matrix(:, 2), proj_plane_corner_matrix(:, 3), 'r', 'FaceAlpha', 0.25, 'DisplayName', 'Plane');   
    
                    % Axes
                    xlabel('a [m]');
                    ylabel('b [m]');
                    zlabel('c [m]');
    
                    axis equal      

                    max_list    = max([ellipsoid_coord_matrix_r(1:num_dim - 1); proj_ellipse_coordinates], [], 1);
                    min_list    = min([ellipsoid_coord_matrix_r(1:num_dim - 1); proj_ellipse_coordinates], [], 1);
                    ampl_list   = max_list - min_list;
    
                    xlim(0.2*ampl_list(1)*[-1, 1] + [min_list(1), max_list(1)]);
                    ylim(0.2*ampl_list(2)*[-1, 1] + [min_list(2), max_list(2)]);
                    
                    view(0, 90);
    
                    set(gca, 'FontSize', 15);
                    set(gca, 'LineWidth', 2);
    
                    hold off    
    
                    % Pause
                    disp('The script will continue and the figure will close upon a button-press');
                    pause();
    
                    close(1);   
                end
            end
end