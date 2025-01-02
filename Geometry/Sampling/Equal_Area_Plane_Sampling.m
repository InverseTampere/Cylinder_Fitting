% The given finite plane is sampled evenly in terms of area.
% The resulting number of points may differ slightly from the given number.

function [plane_sample_matrix, number_samples] = Equal_Area_Plane_Sampling(Plane_Geometry, number_samples, Plot)
    
    %% Inputs %%
        % Plane geometry
        plane_vector_basis  = Plane_Geometry.plane_vector_basis;
        plane_centre        = Plane_Geometry.plane_centre;
        plane_extent_list   = Plane_Geometry.plane_extent_list;

    %% Samples %%
        % 2D samples on the plane
        number_samples_a = sqrt(number_samples * plane_extent_list(1)/plane_extent_list(2));
        number_samples_a = ceil(number_samples_a);
        number_samples_b = sqrt(number_samples * plane_extent_list(2)/plane_extent_list(1));
        number_samples_b = ceil(number_samples_b);
        
        a_steps_list = plane_extent_list(1)/2 * linspace(-1, 1, number_samples_a);
        b_steps_list = plane_extent_list(2)/2 * linspace(-1, 1, number_samples_b);
        
        [a_steps_matrix, b_steps_matrix]    = meshgrid(a_steps_list, b_steps_list); 
        a_steps_list                        = a_steps_matrix(:);
        b_steps_list                        = b_steps_matrix(:);
        number_samples                      = length(a_steps_list);

        proj_plane_sample_matrix            = [a_steps_list, b_steps_list, zeros(number_samples, 1)];

        % Rotated by the plane vector basis
        plane_sample_matrix_t   = (plane_vector_basis' * proj_plane_sample_matrix')';

        % Translated by the plane centre
        plane_sample_matrix     = plane_sample_matrix_t + plane_centre;

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     
            
            hold on
            grid on
            
            % The samples
            scatter3(plane_sample_matrix(:, 1), plane_sample_matrix(:, 2), plane_sample_matrix(:, 3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Samples');

            % The plane
            corner_vector_factors   = {[-1; -1], [-1; 1], [1; 1], [1; -1]};                                                                                     % The planar vectors can either be taken in the positive or negative direction
            Plane_Corner_fun        = @(vector_factors) plane_centre + sum(1/2*plane_extent_list' .* vector_factors .* plane_vector_basis(1 : 2, :), 1);        % Note that the extents were defined as double-sided
            plane_corner_cell       = cellfun(Plane_Corner_fun, corner_vector_factors, 'UniformOutput', false);
            plane_corner_matrix     = vertcat(plane_corner_cell{:});


            patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'b', 'FaceAlpha', 0.25, 'DisplayName', 'Plane');   

             % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            axis equal

            % Viewing angle
            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Pause message
            disp('The samples have been generated. The figure will close and script will end upon a key-press.');
            pause();

            close(1);
        end
end