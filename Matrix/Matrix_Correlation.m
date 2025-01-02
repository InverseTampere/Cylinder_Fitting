% The correlation coefficient of the m*n matrix is computed between the n-variables

function correlation_coeff_matrix = Matrix_Correlation(data_matrix)

    %% Correlation %%
        % The correlation coefficient matrix is initiated with ones, as this is the autocorrelation value
        number_variables            = size(data_matrix, 2);
        correlation_coeff_matrix    = ones(number_variables);
        
        for i = 1 : number_variables
            data_list_i = data_matrix(:, i);
            
            for j = i + 1 : number_variables
                data_list_j = data_matrix(:, j);
                
                % Correlation between variables i and j
                rho_matrix  = corrcoef(data_list_i, data_list_j);
                rho         = rho_matrix(1, 2);                     % The off-diagonal entry is the correlation between i and j
                
                correlation_coeff_matrix(i, j) = rho;
                correlation_coeff_matrix(j, i) = rho;
            end
        end

end