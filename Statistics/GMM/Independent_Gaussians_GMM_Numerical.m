% The 1D Gaussians of each variable within the Gaussian mixture model are computed numerically over the variable data used to create the GMM
% The confidence interval is expected to be given in percent

% Note: The right confidence intervals are given in the first row, the left in the second

function [Variable_Distributions, variable_confidence_interval_matrix, variable_mu_list, variable_STD_list] = Independent_Gaussians_GMM_Numerical(Variable_Data, Confidence_interval)

    %% Independent Gaussians for each variable %%
        % Number of standard deviations corresponding to the confidence interval
        m_STD = sqrt(chi2inv(Confidence_interval / 100, 1));

        % Structure that will contain the data using the same variable names as the data structure
        Variable_Distributions  = struct();

        variable_names          = fieldnames(Variable_Data);
        number_variables        = length(variable_names);

        % The confidence intervals, expected values and standard deviations
        variable_confidence_interval_cell   = cell(1, number_variables);
        variable_mu_cell                    = cell(1, number_variables);
        variable_STD_cell                   = cell(1, number_variables);

        for v = 1 : number_variables
            % This variable's data
            variable    = variable_names{v};
            data_cell   = {Variable_Data.(variable)};
            data_matrix = vertcat(data_cell{:});

            % The expected value, standard deviation and confidence interval
            mu_list             = mean(data_matrix, 1);
            sigma_list          = std(data_matrix, 0, 1);
            confidence_matrix   = mu_list + m_STD * [-1; 1] .* sigma_list;

            % Appended to the structure and cell arrays
            Variable_Distributions.(variable)       = struct('mu', mu_list, 'sigma', sigma_list, 'confidence_interval', confidence_matrix);

            variable_confidence_interval_cell{v}    = confidence_matrix;
            variable_mu_cell{v}                     = mu_list;
            variable_STD_cell{v}                    = sigma_list;
        end

        % Conversion from cell arrays to regular arrays
        variable_confidence_interval_matrix = horzcat(variable_confidence_interval_cell{:});
        variable_mu_list                    = horzcat(variable_mu_cell{:});
        variable_STD_list                   = horzcat(variable_STD_cell{:});

end

