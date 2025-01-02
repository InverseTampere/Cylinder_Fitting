% The correlation coefficients are determined between the fields of the given structure
% Within the structure, the correlation between parameters A and B are saved as Correlation.A.B and Correlation.B.A

function Correlation = Structure_Data_Correlation(Data_Structure, Confidence_interval, Print)
    
    %% Correlation coefficients %%
        % The correlation data is saved in a structure
        variable_names      = fieldnames(Data_Structure);
        number_variables    = length(variable_names);
        corr_struct_start   = [variable_names'; cell(1, number_variables)];
        Correlation         = struct(corr_struct_start{:});
        
        % Statistical significance is determined according to the given confidence interval
        P_threshold = Confidence_interval / 100;
        
        % As the data may be multi-dimensional, correlation is computed between each of the dimensions
        for i = 1 : number_variables
            variable_i      = variable_names{i};
            data_matrix_i   = Data_Structure.(variable_i);
            num_dim_i       = size(data_matrix_i, 2);
            
            for j = i : number_variables
                variable_j      = variable_names{j};
                data_matrix_j   = Data_Structure.(variable_j);
                num_dim_j       = size(data_matrix_j, 2);
                
                % Correlation between the dimensions of the variables' data
                correlation_coeff_matrix    = NaN(num_dim_i, num_dim_j); 
                significance_matrix         = false(num_dim_i, num_dim_j);
                
                for n = 1 : num_dim_i
                    data_list_n = data_matrix_i(:, n);
                    
                    for m = 1 : num_dim_j
                        data_list_m = data_matrix_j(:, m);
                        
                        [rho_matrix, P_matrix] = corrcoef(data_list_n, data_list_m);
                        
                        correlation_coefficient         = rho_matrix(1, 2);                 % Note that the cross-values are off-diagonal
                        correlation_coeff_matrix(n, m)  = correlation_coefficient;
                        
                        P                               = P_matrix(1, 2);
                        statistical_significance        = P > P_threshold;
                        significance_matrix(n, m)       = statistical_significance;
                        
                        % If desired, the results are printed
                        if Print == true
                            if statistical_significance == true
                                fprintf('%s - %s: rho = %.3g, P = %.3g (significant) \n', variable_i, variable_j, correlation_coefficient, P);
                            else
                                fprintf('%s - %s: rho = %.3g, P = %.3g (NOT significant) \n', variable_i, variable_j, correlation_coefficient, P);                                
                            end
                        end
                    end
                end
                
                % Appended to the structure
                Correlation.(variable_i).(variable_j) = struct('correlation_coeff_matrix', correlation_coeff_matrix, 'significance_matrix', significance_matrix);
                Correlation.(variable_j).(variable_i) = struct('correlation_coeff_matrix', correlation_coeff_matrix', 'significance_matrix', significance_matrix');     % Note that the transpose matrices are used
            end
        end

end