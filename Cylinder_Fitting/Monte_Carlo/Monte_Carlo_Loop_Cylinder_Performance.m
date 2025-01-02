% The relative error (in percent) for each individual fitted cylinder within the Monte Carlo loop is determined 
% The delta is negative if the optimiser has a lower error than least-squares

function [LS_Geometry, Fuzzy_Geometry, GMM_Geometry, Delta_Errors, objective_value_list, number_steps_list, number_UU_iter_list] = Monte_Carlo_Loop_Cylinder_Performance(Monte_Carlo_Loop, True_Cylinder_Geometry, LS_Geometry, Fuzzy_Geometry, GMM_Geometry, variable_fields, normalisation_fields)

    %% General performance data %%
        % Objective values at the final step
        Optimiser_Steps_cell    = {Monte_Carlo_Loop.Optimiser_Diagnostics.Optimiser_Steps};
        Objective_Value_fun     = @(Optimiser_Steps) Optimiser_Steps.Objective_value_steps(end);
        objective_value_list    = cellfun(Objective_Value_fun, Optimiser_Steps_cell);

        % Number of steps for each iteration
        Number_Steps_fun        = @(Optimiser_Steps) Optimiser_Steps.num_optimiser_steps;
        number_steps_list       = cellfun(Number_Steps_fun, Optimiser_Steps_cell);

        % Number of uncertainty update iterations for each iteration
        number_UU_iter_cell     = {Monte_Carlo_Loop.Optimiser_Diagnostics.number_UU_iterations};
        number_UU_iter_list     = vertcat(number_UU_iter_cell{:});

    %% Relative errors and uncertainty added to the geometry structures %%    
        % The sets of geometry
        Geometry_sets_cell  = {Fuzzy_Geometry, LS_Geometry, GMM_Geometry};
        geom_set_names      = {'Fuzzy', 'LS', 'GMM'};

        number_geom_sets    = length(geom_set_names);
        number_variables    = length(variable_fields);

        for i = 1 : number_geom_sets
            % Geometry data
            Geometry_Data = Geometry_sets_cell{i};

            % For each variable
            for v = 1 : number_variables
                % The variables' data
                variable        = variable_fields{v};
                variable_mu     = Geometry_Data.(variable).mu;
                variable_sigma  = Geometry_Data.(variable).sigma;
                true_mu         = True_Cylinder_Geometry.(variable);

                % The normalising value
                normalising_field = normalisation_fields{v};

                if isnumeric(normalising_field)                 % Normalised by the value
                    normalising_value = normalising_field;
                else                                            % Normalised by the true geometry
                    normalising_value = abs(True_Cylinder_Geometry.(normalising_field));
                end

                % The relative error
                relative_error                          = (variable_mu - true_mu) / normalising_value * 100;         
                Geometry_Data.(variable).relative_error = relative_error;

                % The relative uncertainty
                relative_sigma                          = variable_sigma / normalising_value * 100;
                Geometry_Data.(variable).relative_sigma = relative_sigma;
            end

            % The new structure is added to the cell
            Geometry_sets_cell{i} = Geometry_Data;
        end

        % The cell array is split
        [Fuzzy_Geometry, LS_Geometry, GMM_Geometry] = Column_Deal(Geometry_sets_cell);

    %% Delta errors %%
        % The difference in error between least squares and the other geometry sets
        Delta_Errors = struct();

        for i = 1 : number_geom_sets
            % Geometry data
            geom_set_name = geom_set_names{i};
            Geometry_Data = Geometry_sets_cell{i};

            if strcmp(geom_set_name, 'LS')          % Least-squares data is skipped
                continue
            end

            % For each variable
            for v = 1 : number_variables
                % The variables' relative error
                variable            = variable_fields{v};
                relative_error      = Geometry_Data.(variable).relative_error;
                LS_relative_error   = LS_Geometry.(variable).relative_error;

                % The delta error                
                delta_error                                 = relative_error - LS_relative_error;
                Delta_Errors.(geom_set_name).(variable)     = delta_error;
            end
        end
end