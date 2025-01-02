% The history of the optimiser is saved, but can have a different length for each Monte Carlo iteration
% In this script, the history is padded s.t. they are equal length, and at each step the mean and standard deviations are computed

function [Fuzzy_Geometry, Least_Squares_Geometry, Objective_Steps, Geometry_Steps, Geometry_Gradient_Steps, Step_Sizes] = Cylinder_Optimisation_History(Monte_Carlo_Loop)
    
    %% Objective value history %%
        % The objective values are converted to a padded matrix
        Optimiser_Steps_cell    = {Monte_Carlo_Loop.Optimiser_Diagnostics.Optimiser_Steps};

        Objective_Steps_fun     = @(Optimiser_Steps) Optimiser_Steps.Objective_value_steps;
        objective_steps_cell    = cellfun(Objective_Steps_fun, Optimiser_Steps_cell, 'UniformOutput', false);

        padding_value           = NaN;            % The matrices are padded with NaN for steps between the length of one loop and the maximum lenght
        objective_steps_matrix  = Padded_Vector_Matrix(objective_steps_cell, padding_value);        
    
        % The number of steps
        Optimiser_Steps_cell    = {Monte_Carlo_Loop.Optimiser_Diagnostics.Optimiser_Steps};
        Number_Steps_fun        = @(Optimiser_Steps) Optimiser_Steps.num_optimiser_steps;
        number_steps_list       = cellfun(Number_Steps_fun, Optimiser_Steps_cell);

        max_num_optim_steps     = max(number_steps_list);
        
        % The mean value and standard deviation
        mean_objective_list     = mean(objective_steps_matrix, 2, 'omitnan');       % Note that NaN values are ignored
        objective_STD_list      = std(objective_steps_matrix, 0, 2, 'omitnan');     % Note that NaN values are ignored
        
        Objective_Steps = struct('mu', mean_objective_list, 'sigma', objective_STD_list, 'matrix', objective_steps_matrix);
        
    %% Optimal and least-squares geometry %%
        % The collection of geometry determined during the Monte Carlo loop
        MC_LS_Geometry      = Monte_Carlo_Loop.MC_Geometry_Data.LS;
        MC_Fuzzy_Geometry   = Monte_Carlo_Loop.MC_Geometry_Data.Optimal;
        
        % The mean values and standard deviations are saved in two structures
        Fuzzy_Geometry          = struct();
        Least_Squares_Geometry  = struct();
        
        geometry_labels         = fieldnames(MC_LS_Geometry);
        num_geo_types           = length(geometry_labels);
        
        for g = 1 : num_geo_types
            % This geometry type's data
            geometry_label      = geometry_labels{g};
            
            geometry_LS_data    = vertcat(MC_LS_Geometry.(geometry_label));
            geometry_opt_data   = vertcat(MC_Fuzzy_Geometry.(geometry_label));
            
            % The mean and standard deviations
            LS_mean     = mean(geometry_LS_data, 1);        % Note that the mean and standard deviations are taken row-wise
            LS_STD      = std(geometry_LS_data, 0, 1);
            
            Opt_mean    = mean(geometry_opt_data, 1);
            Opt_STD     = std(geometry_opt_data, 0, 1);

            % If it is the cylinder direction, the mean is normalised to ensure it is of unit norm length
            if strcmp(geometry_label, 'Cylinder_direction')
                LS_mean     = LS_mean / norm(LS_mean);
                Opt_mean    = Opt_mean / norm(Opt_mean);
            end
            
            % Appended to the two structures
            Least_Squares_Geometry.(geometry_label) = struct('mu', LS_mean, 'sigma', LS_STD);
            Fuzzy_Geometry.(geometry_label)       = struct('mu', Opt_mean, 'sigma', Opt_STD);
        end
        
    %% Optimiser history %%
        % The steps taken by the optimiser in terms of the geometry and its gradients
        Geometry_Steps_fun      = @(Optimiser_Steps) Optimiser_Steps.Geometry_Steps;
        MC_Geometry_Steps_cell  = cellfun(Geometry_Steps_fun, Optimiser_Steps_cell, 'UniformOutput', false);
        MC_Geometry_Steps       = vertcat(MC_Geometry_Steps_cell{:});

        Geo_Grad_Steps_fun      = @(Optimiser_Steps) Optimiser_Steps.Gradient_Steps;
        MC_Geo_Grad_Steps_cell  = cellfun(Geo_Grad_Steps_fun, Optimiser_Steps_cell, 'UniformOutput', false);
        MC_Geo_Grad_Steps       = vertcat(MC_Geo_Grad_Steps_cell{:});

        data_labels     = fieldnames(MC_Geometry_Steps);
        num_data_types  = length(data_labels);
        
        % The mean and standard deviation of the steps are saved per geometry type (gradient) as well as the dimensionality of this data
        Geometry_Steps          = [];
        Geometry_Gradient_Steps = [];
        num_data_dimensions     = 0;
        
        for g = 1 : num_data_types
            % This geometry type's data
            geometry_label                  = data_labels{g};
            geometry_steps_cell             = {MC_Geometry_Steps.(geometry_label)};
            geometry_gradient_steps_cell    = {MC_Geo_Grad_Steps.(geometry_label)};

            % The geometry steps may be multidimensional (i.e. location or vector)
            num_dim             = size(geometry_steps_cell{1}, 2);
            num_data_dimensions = num_data_dimensions + num_dim;
            
            geometry_mean_matrix  = zeros(max_num_optim_steps, num_dim);
            geometry_STD_matrix   = zeros(max_num_optim_steps, num_dim);

            geometry_gradient_mean_matrix   = zeros(max_num_optim_steps, num_dim);
            geometry_gradient_STD_matrix    = zeros(max_num_optim_steps, num_dim);

            for d = 1 : num_dim                
                % This dimension of the cell array is selected
                Dimension_fun                   = @(matrix) matrix(:, d);
                geometry_dim_steps_cell         = cellfun(Dimension_fun, geometry_steps_cell, 'UniformOutput', false);
                geom_gradient_dim_steps_cell    = cellfun(Dimension_fun, geometry_gradient_steps_cell, 'UniformOutput', false);
                
                % Converted to padded matrices
                geometry_steps_matrix           = Padded_Vector_Matrix(geometry_dim_steps_cell, padding_value);
                geometry_gradient_steps_matrix  = Padded_Vector_Matrix(geom_gradient_dim_steps_cell, padding_value);

                % The standard deviation and mean
                geometry_mean_matrix(:, d)      = mean(geometry_steps_matrix, 2, 'omitnan');        % Note that NaN values are ignored
                geometry_STD_matrix(:, d)       = std(geometry_steps_matrix, 0, 2, 'omitnan');      

                geometry_gradient_mean_matrix(:, d) = mean(geometry_gradient_steps_matrix, 2, 'omitnan');
                geometry_gradient_STD_matrix(:, d)  = std(geometry_gradient_steps_matrix, 0, 2, 'omitnan');
            end
            
            % Appended to the structures
            Geometry_Steps.(geometry_label)             = struct('mu', geometry_mean_matrix, 'sigma', geometry_STD_matrix, 'matrix', geometry_steps_matrix);
            Geometry_Gradient_Steps.(geometry_label)    = struct('mu', geometry_gradient_mean_matrix, 'sigma', geometry_gradient_STD_matrix, 'matrix', geometry_gradient_steps_matrix);
        end

        % The step sizes
        Step_Size_fun           = @(Optimiser_Steps) Optimiser_Steps.Step_sizes;
        step_size_cell          = cellfun(Step_Size_fun, Optimiser_Steps_cell, 'UniformOutput', false);
        step_size_matrix        = Padded_Vector_Matrix(step_size_cell, padding_value);

        step_size_mean_list     = mean(step_size_matrix, 2, 'omitnan');         % Note that NaN values are ignored
        step_size_STD_list      = std(step_size_matrix, 0, 2, 'omitnan');

        Step_Sizes              = struct('mu', step_size_mean_list, 'sigma', step_size_STD_list, 'matrix', step_size_matrix);
end