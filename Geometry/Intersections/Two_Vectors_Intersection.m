% The intersection is computed between a given set of 2D vectors and another (constant) 2D vector
% The distances omega and delta are along the set of vectors and constant vector respectively

function [intersection_matrix, omega_list, delta_list] = Two_Vectors_Intersection(start_point_matrix, vector_matrix, constant_start_point, constant_vector, Plot)

    % Normalised vectors
    vector_matrix   = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));
    constant_vector = constant_vector / norm(constant_vector);
    
    % Distances at which intersections occur
    if constant_vector(1) == 0
        omega_list = (constant_start_point(1) - start_point_matrix(:, 1)) ./ vector_matrix(:, 1);
        delta_list = start_point_matrix(:, 2) - constant_start_point(2) + omega_list .* vector_matrix(:, 2);   
    elseif constant_vector(2) == 0
        omega_list = (constant_start_point(2) - start_point_matrix(:, 2)) ./ vector_matrix(:, 2);
        delta_list = start_point_matrix(:, 1) - constant_start_point(1) + omega_list .* vector_matrix(:, 1);
    else        
        delta_list = 1./(1 - constant_vector(2)/constant_vector(1) * vector_matrix(:, 1) ./ vector_matrix(:, 2)) .* ((start_point_matrix(:, 1) - constant_start_point(1))./constant_vector(1) + vector_matrix(:, 1)./(constant_vector(1) * vector_matrix(:, 2)) .* (constant_start_point(2) - start_point_matrix(:, 2)));
        omega_list = (constant_start_point(2) - start_point_matrix(:, 2) + delta_list.*constant_vector(2)) ./ vector_matrix(:, 2);
        
        horizontal_vectors              = vector_matrix(:, 2) == 0;
        delta_list(horizontal_vectors)  = (start_point_matrix(horizontal_vectors, 2) - constant_start_point(2)) / constant_vector(2);
        omega_list(horizontal_vectors)  = constant_start_point(1) - start_point_matrix(horizontal_vectors, 1) + delta_list(horizontal_vectors) * constant_vector(1);
    end

    intersection_matrix = constant_start_point + delta_list .* constant_vector;
    
    % Plot
    if Plot == true
        % Vector colour map
        number_vectors  = size(vector_matrix, 1);
        vector_cmap     = cbrewer('qual', 'Set1', max(number_vectors, 3));
        
        figure(1)
        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
        set(gcf, 'color', [1, 1, 1])     

        hold on
        grid on
        
        % Constant vector
        vector_scaling = 1.2*max(delta_list);
        scatter(constant_start_point(1), constant_start_point(2), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b', 'DisplayName', 'Constant start point');        
        quiver(constant_start_point(1), constant_start_point(2), vector_scaling*constant_vector(1), vector_scaling*constant_vector(2), 'LineWidth', 1, 'color', 'b', 'DisplayName', 'Constant vector');

        for v = 1 : number_vectors
            % This vector's data
            vector_colour   = vector_cmap(v, :);
            start_point     = start_point_matrix(v, :);
            vector          = vector_matrix(v, :);
            
            omega           = omega_list(v);
            intersection    = intersection_matrix(v, :);
            
            scatter(start_point(1), start_point(2), 'MarkerEdgeColor', vector_colour, 'MarkerFaceColor', 'none', 'DisplayName', sprintf('Start %g', v));
            scatter(intersection(1), intersection(2), 'filled', 'MarkerFaceColor', vector_colour, 'DisplayName', sprintf('Intersection %g', v));
            
            vector_scaling  = 1.2*omega; 
            quiver(start_point(1), start_point(2), vector_scaling*vector(1), vector_scaling*vector(2), 'LineWidth', 1, 'color', vector_colour, 'DisplayName', sprintf('Vector %i', v));
        end
        
        % Axes
        axis equal
        xlabel('x [m]');
        ylabel('y [m]');

        % Legend
        legend('show', 'location', 'eastoutside');

        set(gca, 'FontSize', 15);
        set(gca, 'LineWidth', 2);

        hold off  
        
        % Pause message
        disp('The vector intersections have been found. The script will finish and figure will close upon a key-press.');
        pause();
        
        close(1);
    end
end

