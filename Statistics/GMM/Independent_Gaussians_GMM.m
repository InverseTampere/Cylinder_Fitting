% The 1D Gaussians of each variable within the Gaussian mixture model are computed by projecting the mixture's n-D ellipsoids to arrive at a single solution
% The confidence interval is expected to be given in percent

% Note: The Variable_Data structure used to create the GMM is optional. If given, the variable names and their dimensionality are used to create the output structure
% Note: The right confidence intervals are given in the first row, the left in the second

function [Variable_Distributions, variable_confidence_interval_matrix, variable_mu_list, variable_STD_list] = Independent_Gaussians_GMM(GM_Model, Variable_Data, Confidence_interval, convergence_threshold, max_iterations, Print, Diagnostics)

    %% Gaussian mixture model properties %%
        % The properties are taken from the structure
        Mu_matrix           = GM_Model.mu;                      % Expected values
        Weights_list        = GM_Model.ComponentProportion;     % The weight of each Gaussian
        number_components   = GM_Model.NumComponents;           % Number of Gaussians in the model
        number_variables    = GM_Model.NumVariables;            % Dimensionality of each Gaussian
        Sigma_matrix        = GM_Model.Sigma;                   % Covariance matrix
        Shared_Covariance   = GM_Model.SharedCovariance;        % Whether or not the components have the same or different covariance matrices
        
    %% Multivariate normal projection %%
        % As the space is simply the variables themselves, the projection matrix equals the identity matrix and projection point the origin
        vector_matrix   = eye(number_variables);
        origin          = zeros(1, number_variables);

        % The covariance and expected value matrices are changed to cell arrays for cellfun
        if Shared_Covariance == true
            Sigma_matrix = repmat(Sigma_matrix, [1, 1, number_components]);        % To simplify the code, it is repeated if it is shared
        end
        
        Sigma_cell  = mat2cell(Sigma_matrix, number_variables, number_variables, ones(1, number_components));
        Sigma_cell  = squeeze(Sigma_cell);
        
        mu_cell     = mat2cell(Mu_matrix, ones(1, number_components), number_variables);
        
        % Projection
        projection_fun = @(mu, covariance_matrix) Multivariate_Normal_Vector_Projection(vector_matrix, origin, mu, covariance_matrix);

        [proj_mu_cell, ~, proj_extent_cell] = cellfun(projection_fun, mu_cell, Sigma_cell, 'UniformOutput', false);
        proj_mu_cell                        = cellfun(@diag, proj_mu_cell, 'UniformOutput', false);                     % As the vector matrix is an identity matrix, there are only diagonal components
        
        proj_mu_matrix      = horzcat(proj_mu_cell{:});
        proj_sigmae_matrix  = horzcat(proj_extent_cell{:}) / 2;                                                         % Divided by two as the extent is double-sided
        
    %% 1D expected values and standard deviations %%
        % The average weighted values of the projected 1D distributions
        weighted_proj_mu_matrix     = Weights_list .* proj_mu_matrix;
        variable_mu_list            = sum(weighted_proj_mu_matrix, 2)';

        weighted_proj_STD_matrix    = Weights_list .* proj_sigmae_matrix;
        variable_STD_list           = sum(weighted_proj_STD_matrix, 2)';

    %% 1D confidence intervals %%
        % For a mixture model, there is no analytical quantile function
        % The confidence interval is thus found iteratively on the cumulative density function, 
        % where the search interval is derived using the number of standard deviations corresponding to the coverage probability
        P       = Confidence_interval / 100;
        m_STD   = sqrt(chi2inv(P, number_variables));
        
        % The confidence interval leads to the following bounds for the cumulative density
        F_LB            = (1 - P)/2;
        F_UB            = 1 - F_LB;
        F_bounds_list   = [F_LB, F_UB];             % Note that this notation is used s.t. the right- and left-sided bounds are on the correct side of the matrix
        number_bounds   = length(F_bounds_list);
        
        % The metrics of each variable are found separately
        variable_confidence_interval_matrix = zeros(number_bounds, number_variables);

        for v = 1 : number_variables
            % Properties of the projected GMM components for this variable
            mu_list     = proj_mu_matrix(v, :);
            sigmae_list = proj_sigmae_matrix(v, :);
            
            % The value closest to each probability bound is found using iterative sampling due to robustness over the following search interval
            search_lower_bound = min(mu_list - m_STD * sigmae_list);
            search_upper_bound = max(mu_list + m_STD * sigmae_list);
            
            for f = 1 : number_bounds
                F_bound = F_bounds_list(f);
                
                Mixture_CDF_Handle      = @(t) Delta_Cumulative_Density_Function(t);
                [variable_bound, ~, ~]  = Sampling_Function_Minimiser(Mixture_CDF_Handle, search_lower_bound, search_upper_bound, convergence_threshold, max_iterations, Print, Diagnostics);
                
                variable_confidence_interval_matrix(f, v) = variable_bound;
            end
            
            % Diagnostics plot
            if Diagnostics == true
                % The Gaussian probability and cumulative density values
                number_samples          = 1e3;
                search_interval_list    = linspace(search_lower_bound, search_upper_bound, number_samples);
                
                Gaussian_PDF_matrix     = zeros(number_components, number_samples);
                Gaussian_CDF_matrix     = zeros(number_components, number_samples);
                
                for c = 1 : number_components
                    mu      = mu_list(c);
                    sigma   = sigmae_list(c);
                    
                    Gaussian_PDF_list           = normpdf(search_interval_list, mu, sigma);
                    Gaussian_PDF_matrix(c, :)   = Gaussian_PDF_list;
                    
                    Gaussian_CDF_list           = normcdf(search_interval_list, mu, sigma);
                    Gaussian_CDF_matrix(c, :)   = Gaussian_CDF_list;
                end
                
                GMM_PDF_list = sum(Weights_list' .* Gaussian_PDF_matrix, 1);
                GMM_CDF_list = sum(Weights_list' .* Gaussian_CDF_matrix, 1);
                
                % The y-axis of the plot is limited by the Gaussian/GMM values
                pdf_upper_limit = max([Gaussian_PDF_matrix; GMM_PDF_list], [], 'all') + 0.1*(max([Gaussian_PDF_matrix; GMM_PDF_list], [], 'all') - min([Gaussian_PDF_matrix; GMM_PDF_list], [], 'all'));
                            
                figure(v)
                % Size and white background
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])   

                %--% Probability density values %--%
                subplot(1, 2, 1)
                hold on
                grid on

                % The individual Gaussian PDF values
                for c = 1 : number_components
                    Gaussian_PDF_list   = Gaussian_PDF_matrix(c, :);
                    pl_gaussian         = plot(search_interval_list, Gaussian_PDF_list, 'color', 'b', 'LineWidth', 2, 'DisplayName', 'Individual Gaussians');
                    
                    if c > 1
                        pl_gaussian.HandleVisibility = 'Off';
                    end
                end

                % The 1D GMM PDF
                plot(search_interval_list, GMM_PDF_list, 'color', 'm', 'LineWidth', 2, 'DisplayName', '1D GMM');
                
                % The confidence interval
                for f = 1 : number_bounds                
                    variable_bound  = variable_confidence_interval_matrix(f, v);
                    pl_var_bound    = plot(variable_bound * [1,1], [0, pdf_upper_limit], 'LineWidth', 2, 'color', 'r', 'LineStyle', '--', 'DisplayName', 'Confidence interval');
                    
                    if f > 1
                        pl_var_bound.HandleVisibility = 'Off';
                    end
                end
                
                % Axes
                xlabel('t');
                ylabel('f(t)');
                
                xlim([search_lower_bound, search_upper_bound]);
                ylim([0, pdf_upper_limit]);

                % Legend
                legend('show', 'location', 'northoutside');
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
                
                %--% Cumulative density values %--%
                subplot(1, 2, 2)
                hold on
                grid on

                % The individual Gaussian CDF values
                for c = 1 : number_components
                    Gaussian_CDF_list   = Gaussian_CDF_matrix(c, :);
                    pl_gaussian         = plot(search_interval_list, Gaussian_CDF_list, 'color', 'b', 'LineWidth', 2, 'DisplayName', 'Individual Gaussians');
                    
                    if c > 1
                        pl_gaussian.HandleVisibility = 'Off';
                    end
                end

                % The 1D GMM CDF
                plot(search_interval_list, GMM_CDF_list, 'color', 'm', 'LineWidth', 2, 'DisplayName', '1D GMM');
                
                % The confidence interval
                for f = 1 : number_bounds                
                    variable_bound  = variable_confidence_interval_matrix(f, v);
                    pl_var_bound    = plot(variable_bound * [1,1], [0, 1], 'LineWidth', 2, 'color', 'r', 'LineStyle', '--', 'DisplayName', 'Confidence interval');
                    
                    F_bound         = F_bounds_list(f);
                    pl_F_bound      = plot([search_lower_bound, search_upper_bound], F_bound  * [1,1], 'LineWidth', 2, 'color', 'r', 'LineStyle', '--', 'DisplayName', 'Confidence interval');
                    pl_F_bound.HandleVisibility = 'Off';
                    
                    if f > 1
                        pl_var_bound.HandleVisibility = 'Off';
                    end
                end
                
                % Axes
                xlabel('t');
                ylabel('F(t)');
                
                xlim([search_lower_bound, search_upper_bound]);
                ylim([0, 1]);

                % Legend
                legend('show', 'location', 'northoutside');
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
                
                % Paused for analysis
                disp('The figure will close and the script will continue when a key is pressed.');
                pause();
                close(v);                
            end
        end
        
    %% Uncertainty structure %%
        % If the input data structure is given, the uncertainty structure can be generated
        if ~isempty(Variable_Data)
            variable_names_cell         = fieldnames(Variable_Data); 
            number_variable_types       = length(variable_names_cell);
            uncertainty_struct_start    = [variable_names_cell'; cell(1, number_variable_types)];            
            Variable_Distributions      = struct(uncertainty_struct_start{:});

            start_ind = 1;     % The variables may be multi-dimensional, so the correct column is selected this way

            for v = 1 : number_variable_types
                % The number of dimensions this variable occupies
                variable_name   = variable_names_cell{v};
                example_data    = Variable_Data(1).(variable_name);         % The first entry is used just as an example
                num_dim         = length(example_data);

                end_ind         = start_ind + num_dim - 1;
                
                % This variable's expected value, standard deviation and confidence interval
                expected_value      = variable_mu_list(start_ind : end_ind);
                standard_deviation  = variable_STD_list(start_ind : end_ind);
                confidence_interval = variable_confidence_interval_matrix(:, start_ind : end_ind);
                
                Variable_Distributions.(variable_name) = struct('mu', expected_value, 'sigma', standard_deviation, 'confidence_interval', confidence_interval);
                
                % The start index is updated
                start_ind = end_ind + 1;
            end            
        end
        
    %% 1-D Gaussian mixture's cumulative density function %%
    function delta_F = Delta_Cumulative_Density_Function(t) 
        % The cumulative density of this point t for each mixture's component
        F_components    = 1/2 * (1 + erf((t - mu_list) ./ (sqrt(2)*sigmae_list)));
        
        % The over-all weighted cumulative density
        F               = sum(Weights_list .* F_components);
        
        % As the minimum of the function is found, the absolute delta w.r.t. the searched for bound is used
        delta_F         = abs(F - F_bound);
    end        
end

