% The over-all weighted expected values and covariance matrix are computed for the given Gaussian mixture model

% The data structure is optional. When given, a structure is returned containing the results for each variable
% Note that in the structure, the covariance between parameters A and B is given as GMM_Properties.A.B and GMM_Properties.B.A

function [Weighted_expected_values_list, Weighted_STD_list, Weighted_covariance_matrix, GMM_Variable_Properties] = Weighted_GMM_Properties(GM_Model, Data_Structure)
    
    %% GMM properties %%
        Mu_matrix               = GM_Model.mu;                       % Expected values
        Sigma_matrix            = GM_Model.Sigma;                    % Covariance matrix
        Shared_Covariance       = GM_Model.SharedCovariance;         % Whether or not the components have the same or different covariance matrices
        Component_weights_list  = GM_Model.ComponentProportion;      % The weight of each Gaussian
        number_components       = GM_Model.NumComponents;            % Number of Gaussians in the model
        number_variables        = GM_Model.NumVariables;             % The number of variables of the Gaussians

    %% The expected values
        % The expected values are simply the weighted sum
        Weighted_expected_values_list = sum(Component_weights_list' .* Mu_matrix, 1);
        
    %% The covariance matrix %%
        % If they all share the same covariance matrix, the final matrix is identical
        if Shared_Covariance == true 
            Weighted_covariance_matrix = Sigma_matrix;
            
        % It has to be computed if they are independent
        else
            Weighted_covariance_matrix = zeros(number_variables);
            
            for i = 1 : number_variables
                for j = 1 : i
                    % The deviation term
                    Deviation_list = zeros(1, number_components);
                    
                    for c = 1 : number_components
                        % The component's properties for these variables
                        weight  = Component_weights_list(c);
                        mu_i    = Mu_matrix(c, i);
                        mu_j    = Mu_matrix(c, j);
                        
                        sigma   = Sigma_matrix(i, j, c);
                        
                        % The deviation
                        deviation           = weight * (sigma + mu_i*mu_j);
                        Deviation_list(c)   = deviation;
                    end
                    
                    Deviation_term  = sum(Deviation_list);
                    
                    % The mean term
                    Mean_term       = sum(Component_weights_list' .* Mu_matrix(:, i)) * sum(Component_weights_list' .* Mu_matrix(:, j));
                    
                    % The covariance
                    Covariance      = Deviation_term - Mean_term;
                    
                    Weighted_covariance_matrix(i, j) = Covariance;
                    Weighted_covariance_matrix(j, i) = Covariance;
                end
            end
        end

        % The standard deviations equal the square root of the diagonal
        Weighted_STD_list = sqrt(diag(Weighted_covariance_matrix))';
        
    %% Structure %%
        if ~isempty(Data_Structure)
            % The expected values and covariance data are saved in a structure
            variable_names          = fieldnames(Data_Structure);
            number_variables        = length(variable_names);
            properties_start        = [variable_names'; cell(1, number_variables)];
            GMM_Variable_Properties = struct(properties_start{:});

            % As the data may be multi-dimensional, covariance is similarly multidimensional and the right entries are selected by counting rows and columns
            row_start = 1;

            for i = 1 : number_variables
                variable_i      = variable_names{i};
                data_matrix_i   = Data_Structure.(variable_i);
                num_dim_i       = size(data_matrix_i, 2);

                row_end         = row_start + num_dim_i - 1;
                column_start    = 1;

                % The expected value and standard deviation
                GMM_Variable_Properties.(variable_i).mu     = Weighted_expected_values_list(row_start : row_end);
                GMM_Variable_Properties.(variable_i).sigma  = Weighted_STD_list(row_start : row_end);

                for j = i : number_variables
                    variable_j      = variable_names{j};
                    data_matrix_j   = Data_Structure.(variable_j);
                    num_dim_j       = size(data_matrix_j, 2);

                    column_end      = column_start + num_dim_j - 1;

                    % Covariance between the dimensions of the variables' data
                    covariance_matrix_ij = Weighted_covariance_matrix(row_start : row_end, column_start : column_end);                

                    % Appended to the structure
                    GMM_Variable_Properties.(variable_i).(variable_j).Sigma = covariance_matrix_ij;
                    GMM_Variable_Properties.(variable_j).(variable_i).Sigma = covariance_matrix_ij';     % Note that the transpose matrix is used

                    % Update the column start point
                    column_start = column_end + 1;
                end

                % Update the row start point
                row_start = row_end + 1;
            end
            
        else
            GMM_Variable_Properties = [];
        end
        
end