% This script creates samples that are equally spaced and covering an n-dimensional hypercuboidal volume

% If the total number of samples is provided instead of the list for each dimension is not provided, they are determined according to the length of the vertices
% Note that the final number of samples usually differs slightly from the original input in that case

function [cuboid_sample_matrix, number_samples, number_samples_dim_list, delta_volume] = Equal_Volume_Hypercuboid_Sampler(cuboid_centre, cuboid_vertex_lengths, cuboid_coordinate_frame, number_samples, number_samples_dim_list)
        
    %% Sampling the hypercuboid %%
        % Volume of the cuboid
        num_dim     = length(cuboid_centre);
        V_cuboid    = prod(cuboid_vertex_lengths);

        % Number of samples in each dimension
        if isempty(number_samples_dim_list)
            number_samples_dim_list = number_samples^(1/num_dim) * cuboid_vertex_lengths / mean(cuboid_vertex_lengths);
            number_samples_dim_list = ceil(number_samples_dim_list);
        end

        % Delta volume per sample and total number of samples
        number_samples  = prod(number_samples_dim_list);
        delta_volume    = V_cuboid / number_samples;
        
        % Structure containing the samples for each dimension
        Samples_Struct = struct('samples', []);

        for d = 1 : num_dim
            % This dimension's vertex length and number of samples
            dim_vertex_length   = cuboid_vertex_lengths(d);
            number_samples_dim  = number_samples_dim_list(d);

            % The generated samples are saved in the structure
            dim_samples_list            = -dim_vertex_length/2 + dim_vertex_length*linspace(0, 1, number_samples_dim);
            Samples_Struct(d).samples   = dim_samples_list;
        end
        
        % Full matrix of resulting samples
        [cuboid_sample_cell{1: num_dim}] = ndgrid(Samples_Struct(:).samples);

        Vertcat_fun             = @(cuboid_samples) cuboid_samples(:);
        cuboid_sample_cell      = cellfun(Vertcat_fun, cuboid_sample_cell, 'UniformOutput', false);

        cuboid_sample_matrix    = horzcat(cuboid_sample_cell{:});

    %% Transformation to original coordinate frame %%
        % First rotated
        cuboid_sample_matrix = (cuboid_coordinate_frame * cuboid_sample_matrix')';

        % Then translated
        cuboid_sample_matrix = cuboid_sample_matrix + cuboid_centre;

end