% This function determines the intersection between an ellipsoid and a plane.
% It can either be an ellipse, point or not exist.

% If it does not exist, NaN entries are returned for all that plane's outputs.
% If it does exist, but is a point, the axes are the identity matrix and radii is 0

function [ellipse_centre_matrix, ellipse_axes_cell, ellipse_radii_matrix] = Plane_Ellipsoid_Intersection(plane_point_matrix, plane_normal_vector_matrix, ellipsoid_centre, ellipsoid_radii, ellipsoid_axes)

    %% Manual inputs %%
        Plot = false;        % [true, false] Produces a plot of the plane, ellipsoid and their intersecting ellipse (if it exists)

    %% Planar vectors %%
        % Conversion of plane normal vectors to cell array
        [num_planes, num_dim]       = size(plane_normal_vector_matrix);
        plane_normal_vector_cell    = mat2cell(plane_normal_vector_matrix, ones(1, num_planes), num_dim);

        % Planar vectors
        Planar_Vector_fun                                   = @(normal_vector) Vector_Based_Rotation(ellipsoid_centre, normal_vector, ellipsoid_centre);        % Note that the ellipsoid centre is chosen arbitrarily
        [~, plane_vector_basis_cell, planar_vectors_cell]   = cellfun(Planar_Vector_fun, plane_normal_vector_cell, 'UniformOutput', false);

    %% Coordinate system transformation %%        
        % The coordinate system is changed s.t. it is centered on the ellipsoid        
        plane_point_matrix_t    = plane_point_matrix - ellipsoid_centre;
        
        % Aligned with the ellipsoid axes and 'squeezed' s.t. the ellipsoid becomes a unit sphere
        squeezing_matrix        = eye(num_dim) ./ ellipsoid_radii;
        Transformation_fun      = @(vector) (squeezing_matrix * ellipsoid_axes * vector')';

        plane_point_matrix_s    = Transformation_fun(plane_point_matrix_t);
        planar_vectors_cell_s   = cellfun(Transformation_fun, planar_vectors_cell, 'UniformOutput', false);

        % The associated normal vectors
        Normal_Vector_fun       = @(planar_vectors) cross(planar_vectors(1, :) / norm(planar_vectors(1, :)), planar_vectors(2, :) / norm(planar_vectors(2, :)));
        normal_vector_cell_s    = cellfun(Normal_Vector_fun, planar_vectors_cell_s, 'UniformOutput', false);
        normal_vector_matrix_s  = vertcat(normal_vector_cell_s{:});
        
        % The squeezed plane's vector bases
        [~, plane_vector_basis_cell_s, ~] = cellfun(Planar_Vector_fun, normal_vector_cell_s, 'UniformOutput', false);

    %% Unit sphere intersections %%                      
        % This means that the intersection is a simple circle
        unit_sphere_centre = zeros(1, num_dim);
        unit_sphere_radius = 1;

        [circle_centre_matrix, circle_radius_list] = Plane_Sphere_Intersection(plane_point_matrix_s, normal_vector_matrix_s, unit_sphere_centre, unit_sphere_radius);
                
    %% Transformation to intersecting ellipses %%
        % The ellipse centre simply follows from a reverse transformation of the circle centre
        Reverse_Transformation_fun  = @(matrix) (ellipsoid_axes' * (squeezing_matrix \ matrix'))' + ellipsoid_centre;  
        ellipse_centre_matrix       = Reverse_Transformation_fun(circle_centre_matrix);

        % The axes and radii are found by least-squares fitting the intersecting ellipse, which requires 5 circle points
        number_samples                              = 5;
        unit_radius                                 = 1;
        [x_unit_circle_list, y_unit_circle_list]    = Circle_Coordinates(unit_radius, 0, 0, number_samples);
        unit_circle_samples                         = [x_unit_circle_list', y_unit_circle_list', zeros(number_samples, 1)];
    
        Circle_Sampling_fun = @(circle_radius, circle_centre, plane_vector_basis_s) (plane_vector_basis_s' * circle_radius*unit_circle_samples')' + circle_centre;
        
        circle_radius_cell  = num2cell(circle_radius_list);
        circle_centre_cell  = mat2cell(circle_centre_matrix, ones(1, num_planes), num_dim);
        circle_samples_cell = cellfun(Circle_Sampling_fun, circle_radius_cell, circle_centre_cell, plane_vector_basis_cell_s, 'UniformOutput', false);
    
        % They are transformed to the original coordinate frame such that they are ellipse samples
        ellipse_samples_cell        = cellfun(Reverse_Transformation_fun, circle_samples_cell, 'UniformOutput', false);
    
        % They are projected onto the plane s.t. they are 2D
        Proj_Diagnostics            = false;
        Matrix_to_Plane_Projection  = @(matrix, plane_vector_basis, plane_point) Point_to_Plane_Projection(matrix, plane_vector_basis, plane_point, Proj_Diagnostics);

        plane_point_cell                    = mat2cell(plane_point_matrix, ones(1, num_planes), num_dim);
        [~, proj_ellipse_samples_cell, ~]   = cellfun(Matrix_to_Plane_Projection, ellipse_samples_cell, plane_vector_basis_cell, plane_point_cell, 'UniformOutput', false);
    
        % The projected intersecting ellipse is fit exactly
        Ellipse_Fitting_fun                             = @(proj_samples_matrix) Ellipse_Least_Squares_Fitting(proj_samples_matrix);       
        [~, ellipse_radii_cell, proj_ellipse_axes_cell] = cellfun(Ellipse_Fitting_fun, proj_ellipse_samples_cell, 'UniformOutput', false); 
        ellipse_radii_matrix                            = vertcat(ellipse_radii_cell{:});
    
        % Transformed back to the original coordinate frame
        Plane_Deprojection  = @(plane_vector_basis, matrix) (plane_vector_basis' * [matrix, zeros(num_dim - 1, 1)]')';
        ellipse_axes_cell   = cellfun(Plane_Deprojection, plane_vector_basis_cell, proj_ellipse_axes_cell, 'UniformOutput', false);
        
    %% Plot %%
        if Plot == true && num_dim == 3     % A plot is only shown if the input data is 3D
            % Number of points into which the geometry is discretised
            num_coord = 1e3;
            
            % Unit sphere coordinates
            [x_unit_sphere, y_unit_sphere, z_unit_sphere] = sphere(num_coord);
            
            % Ellipsoid coordinates            
            [ellipsoid_coord_matrix, num_coord] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, num_coord);

            % Colour map for the planes            
            plane_cmap  = cbrewer('qual', 'Set1', max(num_planes, 3));
            plane_cmap  = max(plane_cmap, 0);
            plane_cmap  = min(plane_cmap, 1);
            
            % Plot
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            %--% Squeezed coordinate frame %--%
            subplot(1, 2, 1)
            hold on
            grid on
            
            % Unit sphere
            surf(x_unit_sphere, y_unit_sphere, z_unit_sphere, 'EdgeColor', 'none', 'FaceColor', 'b', 'Facealpha', 0.5, 'DisplayName', 'Unit sphere');
            
            % Planes and intersecting circles
            for p = 1 : num_planes
                colour = plane_cmap(p, :);
                
                % Circle
                circle_radius   = circle_radius_list(p);
                circle_centre   = circle_centre_matrix(p, :);
                
                normal_vector_s             = normal_vector_cell_s{p};
                [~, ~, planar_vectors_s]    = Vector_Based_Rotation(circle_centre, normal_vector_s, circle_centre);
                
                theta_list          = linspace(0, 2*pi, num_coord)';
                circle_coord_matrix = circle_centre + circle_radius * (cos(theta_list) .* planar_vectors_s(1, :) + sin(theta_list) .* planar_vectors_s(2, :));
                pl_ci               = plot3(circle_coord_matrix(:, 1), circle_coord_matrix(:, 2), circle_coord_matrix(:, 3), 'LineWidth', 2, 'color', colour, 'DisplayName', 'Intersecting circle');

                % Circle samples
                circle_samples_matrix   = circle_samples_cell{p};
                sc_cs                   = scatter3(circle_samples_matrix(:, 1), circle_samples_matrix(:, 2), circle_samples_matrix(:, 3), 'MarkerFaceColor', colour, 'MarkerEdgeColor', 'k', 'DisplayName', 'Circle samples');

                % Plane
                normal_vector   = normal_vector_matrix_s(p, :);
                plane_point     = plane_point_matrix_s(p, :);
                
                sphere_data_matrix  = repmat([-1; 1], [1, 3]);
                plane_corner_matrix = Plane_Corner_Points(normal_vector, plane_point, [plane_point; sphere_data_matrix]);
                pa_pl = patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), colour, 'FaceAlpha', 0.5, 'DisplayName', 'Plane');   

                if p > 1
                    pl_ci.HandleVisibility = 'Off';
                    sc_cs.HandleVisibility = 'Off';
                    pa_pl.HandleVisibility = 'Off';
                end                
            end
            
            % Axes
            xlabel('u [-]');
            ylabel('v [-]');
            zlabel('w [-]');

            xlim([-3/2, 3/2]);
            ylim([-3/2, 3/2]);
            zlim([-3/2, 3/2]);

            axis equal    

            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
            
            %--% Original coordinate frame %--%
            subplot(1, 2, 2)
            hold on
            grid on

            % Ellipsoid
            x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(num_coord) * [1, 1]);
            y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(num_coord) * [1, 1]);
            z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(num_coord) * [1, 1]);

            surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Ellipsoid');

            % The ellipsoid axes            
            for e = 1 : num_dim
                ellipsoid_axis = ellipsoid_radii(e) * ellipsoid_axes(e, :);   

                e_ax = plot3(ellipsoid_centre(1) + [0, ellipsoid_axis(1)], ellipsoid_centre(2) + [0, ellipsoid_axis(2)], ellipsoid_centre(3) + [0, ellipsoid_axis(3)], 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Ellipsoid axes');

                if e > 1
                    e_ax.HandleVisibility = 'Off';
                end
            end

            % Planes and intersecting ellipses
            for p = 1 : num_planes
                colour = plane_cmap(p, :);
                
                % Ellipse
                ellipse_centre          = ellipse_centre_matrix(p, :);
                ellipse_axes            = ellipse_axes_cell{p};
                ellipse_radii           = ellipse_radii_matrix(p, :);
                ellipse_coord_matrix    = Ellipse_Coordinate_Generator(ellipse_centre, ellipse_axes, ellipse_radii, num_coord);
    
                pl_el = plot3(ellipse_coord_matrix(:, 1), ellipse_coord_matrix(:, 2), ellipse_coord_matrix(:, 3), 'LineWidth', 2, 'color', colour, 'DisplayName', 'Intersecting ellipse');

                % Its samples
                ellipse_samples_matrix  = ellipse_samples_cell{p};
                sc_es                   = scatter3(ellipse_samples_matrix(:, 1), ellipse_samples_matrix(:, 2), ellipse_samples_matrix(:, 3), 'MarkerFaceColor', colour, 'MarkerEdgeColor', 'k', 'DisplayName', 'Ellipse samples');

                % Axes
                for v = 1 : 2
                    ellipse_radius  = ellipse_radii_matrix(p, v); 
                    ellipse_axis    = ellipse_radius * ellipse_axes(v, :);
                    pl_elax = plot3(ellipse_centre(1) + [0, ellipse_axis(1)], ellipse_centre(2) + [0, ellipse_axis(2)], ellipse_centre(3) + [0, ellipse_axis(3)], 'LineWidth', 2, 'color', colour, 'DisplayName', 'Int. ellipse axis');
                    
                    if v > 1 || p > 1
                        pl_elax.HandleVisibility = 'Off';
                    end
                end
                    
                % Plane
                normal_vector   = plane_normal_vector_matrix(p, :);
                plane_point     = plane_point_matrix(p, :);
                
                plane_corner_matrix = Plane_Corner_Points(normal_vector, plane_point, [ellipsoid_coord_matrix; plane_point]);
                pa_pl = patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), colour, 'FaceAlpha', 0.5, 'DisplayName', 'Plane');   
                
                if p > 1
                    pl_el.HandleVisibility = 'Off';
                    sc_es.HandleVisibility = 'Off';
                    pa_pl.HandleVisibility = 'Off';
                end
            end

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            max_list    = max(ellipsoid_coord_matrix, [], 1);
            min_list    = min(ellipsoid_coord_matrix, [], 1);
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
            
            % Time to look at the figure
            disp('The figure is closed when a button is pressed, and the script continues.');
            pause();
            close(1);
        end
        
end