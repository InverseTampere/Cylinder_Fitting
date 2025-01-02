% The correlation coefficients between each of the variables of the Gaussian mixture model can be computed from the over-all covariance matrix
% The data structure is an optional input. When given the correlation coefficients between variables A and B are saved as Correlation.A.B and Correlation.B.A

function [correlation_coeff_matrix, Correlation] = Correlation_GMM(GM_Model, Data_Structure)
    
    %% Correlation coefficients %%
        % The over-all covariance matrix of the Gaussian mixture model
        [~, ~, covariance_matrix, ~] = Weighted_GMM_Properties(GM_Model, []);
    
        % The correlation coefficient between variables i and j are then Sigma_ij / sqrt(Sigma_ii * Sigma_jj)
        variance_list               = diag(covariance_matrix);
        correlation_coeff_matrix    = covariance_matrix ./ sqrt(variance_list' .* variance_list);
        
    %% Structure %%
        % The correlation data is saved in a structure
        variable_names      = fieldnames(Data_Structure);
        number_variables    = length(variable_names);
        corr_struct_start   = [variable_names'; cell(1, number_variables)];
        Correlation         = struct(corr_struct_start{:});
        
        % As the data may be multi-dimensional, correlation is similarly multidimensional and the right entries are selected by counting rows and columns
        row_start = 1;
        
        for i = 1 : number_variables
            variable_i      = variable_names{i};
            data_matrix_i   = Data_Structure.(variable_i);
            num_dim_i       = size(data_matrix_i, 2);
            
            row_end         = row_start + num_dim_i - 1;
            column_start    = 1;
            
            for j = i : number_variables
                variable_j      = variable_names{j};
                data_matrix_j   = Data_Structure.(variable_j);
                num_dim_j       = size(data_matrix_j, 2);
                
                column_end      = column_start + num_dim_j - 1;
                
                % Correlation between the dimensions of the variables' data
                correlation_coeff_matrix_ij = correlation_coeff_matrix(row_start : row_end, column_start : column_end);                
                
                % Appended to the structure
                Correlation.(variable_i).(variable_j) = struct('correlation_coeff_matrix', correlation_coeff_matrix_ij);
                Correlation.(variable_j).(variable_i) = struct('correlation_coeff_matrix', correlation_coeff_matrix_ij');     % Note that the transpose matrix is used
                
                column_start    = column_end + 1;
            end
            
            row_start = row_end + 1;
        end
end