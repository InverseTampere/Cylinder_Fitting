% This script creates samples that are equally spaced and covering an elliptic cylinder
% The ellipse is described by its two axes and radii. It is elongated by the cylinder length, normal to the cross-section

% Note that the number of samples for each dimension is roughly the cubic root of the given number of samples
% The total number of samples may differ slightly from the initially given value

function [sample_coord_matrix, number_samples, number_axis_samples, number_ellipse_samples] = Equal_Volume_Elliptic_Cylinder_Sampler(ellipse_radii, ellipse_axes, el_cyl_length, el_cyl_centroid, number_axis_samples, number_ellipse_samples)
        
    %% Circular samples %%
        % The distance along each of the ellipse axes for equal area sampling
        circle_radius   = 1;
        circle_centre   = zeros(1, 3);
        z_vector        = [0, 0, 1];
        
        [ellipse_norms_matrix, number_ellipse_samples] = Equal_Area_Circular_Sampler(circle_radius, circle_centre, z_vector, number_ellipse_samples);   % Note that the number of elliptical samples may change

        % The resulting samples
        ellipse_samples_matrix = el_cyl_centroid + ellipse_radii(1) * ellipse_norms_matrix(:, 1) .* ellipse_axes(1, :) + ellipse_radii(2) * ellipse_norms_matrix(:, 2) .* ellipse_axes(2, :);     
            
    %% Extension along cylinder axis %%        
        % To ensure the centre is included, the number of axis samples is odd
        if mod(number_axis_samples, 2) == 0
            number_axis_samples = number_axis_samples + 1;
        end
    
        % The height at each step
        height_list = el_cyl_length/2 * linspace(-1, 1, number_axis_samples);
               
        % The vector at each of these heights is added to the ellipse coordinates
        axis_vector          = cross(ellipse_axes(1, :), ellipse_axes(2, :));
        axis_vector_matrix   = height_list' * axis_vector / norm(axis_vector);
        
        height_adjuster     = @(axis_vector) {ellipse_samples_matrix + axis_vector'};
        sample_coord_cell   = splitapply(height_adjuster, axis_vector_matrix', 1 : number_axis_samples);
        sample_coord_matrix = vertcat(sample_coord_cell{:});
        
        % The total number of samples
        number_samples      = number_axis_samples * number_ellipse_samples;

end