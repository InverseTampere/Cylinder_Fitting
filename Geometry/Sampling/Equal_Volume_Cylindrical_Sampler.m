% This script creates samples that are equally spaced and covering a cylindrical volume
% The cylinder is described by its radius, length, vector and centroid

% Note that the number of samples for each dimension is roughly the cubic root of the given number of samples
% The total number of samples may differ slightly from the initially given value

function [sample_coord_matrix, number_samples, number_axis_samples, number_circle_samples] = Equal_Volume_Cylindrical_Sampler(cyl_radius, cyl_length, cyl_vector, cyl_centroid, number_samples)
        
    %% Circular samples %%
        % The number of samples in each circle is the cubic root squared
        number_circle_samples = number_samples^(2/3);
    
        [circle_coord_matrix, number_circle_samples] = Equal_Area_Circular_Sampler(cyl_radius, cyl_centroid, cyl_vector, number_circle_samples);
        
    %% Extension along cylinder axis %%
        % The number of samples along the cylinder axis and resulting total number of samples
        number_axis_samples = round(number_samples / number_circle_samples);
        number_samples      = number_axis_samples * number_circle_samples;
        
        % The resulting height at each step
        height_list = cyl_length/2 * linspace(-1, 1, number_axis_samples);
               
        % The vector at each of these heights is added to the circle coordinates
        cyl_vector          = cyl_vector / norm(cyl_vector);
        cyl_vector_matrix   = height_list' * cyl_vector;
        
        height_adjuster     = @(cyl_vector) {circle_coord_matrix + cyl_vector'};
        sample_coord_cell   = splitapply(height_adjuster, cyl_vector_matrix', 1 : number_axis_samples);
        sample_coord_matrix = vertcat(sample_coord_cell{:});

end