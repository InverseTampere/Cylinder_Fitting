% The confidence intervals of the Gaussian mixture model variables are given,
% taking the outer-most extrema of the projected n-D ellipsoid to arrive at a single solution
% This ensures that the confidence interval is always <= 1 - alpha
% The standard deviation is the average for the left- and right-sided values

% The Variable_Data structure used to create the GMM is optional. If given, the variable names and their dimensionality are used to create the Variable_Uncertainty structure

function [Variable_Uncertainty, variable_bounds_matrix, variable_mu_list, variable_STD_list] = Confidence_Intervals_GMM(GM_Model, Variable_Data, Confidence_interval, convergence_threshold, max_iterations, Print)

    %% Gaussian mixture model properties %%
        % The properties are taken from the structure
        Mu_matrix           = GM_Model.mu;                      % Expected values
        Weights_list        = GM_Model.ComponentProportion;     % The weight of each Gaussian
        number_components   = GM_Model.NumComponents;           % Number of Gaussians in the model
        number_variables    = GM_Model.NumVariables;            % Dimensionality of each Gaussian
        Sigma_matrix        = GM_Model.Sigma;                   % Covariance matrix
        Shared_Covariance   = GM_Model.SharedCovariance;        % Whether or not the components have the same or different covariance matrices
        
    %% Corresponding ellipsoids %%
        % The covariance and expected value matrices are changed to cell arrays
        if Shared_Covariance == true
            Sigma_matrix    = repmat(Sigma_matrix, [1, 1, number_components]);        % To simplify the code, it is repeated if it is shared
        end
        
        Sigma_cell  = mat2cell(Sigma_matrix, number_variables, number_variables, ones(1, number_components));
        Sigma_cell  = squeeze(Sigma_cell);
        
        ellipsoid_centre_cell = mat2cell(Mu_matrix, ones(1, number_components), number_variables);
        
        % The ellipsoid axes and radii are determined from the covariance matrices
        [ellipsoid_axes_cell, eigenvalue_cell]  = cellfun(@eig, Sigma_cell, 'UniformOutput', false);
        
        ellipsoid_radii_fun     = @(eigenvalue_matrix) sqrt(diag(eigenvalue_matrix))';                      % Note that the eigenvalues are only on the diagonal
        ellipsoid_sigmae_cell   = cellfun(ellipsoid_radii_fun, eigenvalue_cell, 'UniformOutput', false);    % And that at this stage, the radii equal the sigmaes
        
    %% Ellipsoid projections %%
        % As the space is simply the variables themselves, the projection matrix equals the identity matrix and projection point the origin
        vector_matrix   = eye(number_variables);
        origin          = zeros(1, number_variables);
    
        % The projected uncertainties and expected values
        projection_fun  = @(ellipsoid_centre, ellipsoid_sigmae, ellipsoid_axes) Ellipsoid_to_Vector_Projection(vector_matrix, origin, ellipsoid_centre, ellipsoid_sigmae, ellipsoid_axes);

        [proj_mu_cell, ~, proj_sigmae_cell] = cellfun(projection_fun, ellipsoid_centre_cell, ellipsoid_sigmae_cell, ellipsoid_axes_cell, 'UniformOutput', false);
        proj_mu_cell                        = cellfun(@diag, proj_mu_cell, 'UniformOutput', false);         % As the vector matrix is an identity matrix, there are only diagonal components
        
        proj_mu_matrix      = horzcat(proj_mu_cell{:});
        proj_sigmae_matrix  = horzcat(proj_sigmae_cell{:});
        
    %% Confidence interval bounds %%
        % For a mixture model, there is no analytical quantile function
        % The confidence interval is thus found using golden-section search on the cumulative density function, 
        % where the search interval is derived using the number of standard deviations corresponding to the coverage probability
        P       = Confidence_interval / 100;
        m_STD   = sqrt(chi2inv(P, number_variables));
        
        % The confidence interval leads to the following bounds for the cumulative density
        F_LB            = (1 - P)/2;
        F_UB            = 1 - F_LB;
        F_bounds_list   = [F_UB; F_LB];             % Note that this notation is used s.t. the upper bound is above the lower bound
        number_bounds   = length(F_bounds_list);
        
        % The metrics of each variable are found separately
        variable_bounds_matrix  = zeros(number_bounds, number_variables);

        for v = 1 : number_variables
            % Properties of the projected GMM components for this variable
            mu_list     = proj_mu_matrix(v, :);
            sigmae_list = proj_sigmae_matrix(v, :);
            
            % The value closest to each probability bound is found using golden-section search over the following search interval
            search_lower_bound = min(mu_list - m_STD * sigmae_list);
            search_upper_bound = max(mu_list + m_STD * sigmae_list);
            
            for f = 1 : number_bounds
                F_bound = F_bounds_list(f);
                
                Mixture_CDF_Handle  = @(t) Delta_Cumulative_Density_Function(t);
                [variable_bound, ~] = Golden_Section_Search(Mixture_CDF_Handle, search_lower_bound, search_upper_bound, convergence_threshold, max_iterations, Print);
                
                variable_bounds_matrix(f, v) = variable_bound;
            end
        end
        
        % Weighted expected value over the projected distributions
        variable_mu_list = sum(Weights_list .* proj_mu_matrix, 2);

        % Single-sided standard deviation
        variable_interval_halfwidths        = abs(variable_bounds_matrix - variable_mu_list);
        variable_single_sided_STD_matrix    = variable_interval_halfwidths / m_STD;
        
        % The average values are used such that a Gaussian is described, as skew is impractically complex to include
        variable_STD_list = sum(variable_single_sided_STD_matrix, 1) / 2;
        
    %% Uncertainty structure %%
        % If the input data structure is given, the uncertainty structure can be generated
        if ~isempty(Variable_Data)
            variable_names_cell         = fieldnames(Variable_Data); 
            uncertainty_struct_start    = [variable_names_cell'; cell(1, number_variables)];            
            Variable_Uncertainty        = struct(uncertainty_struct_start{:});

            start_ind = 1;     % The variables may be multi-dimensional, so the correct column is selected this way

            for v = 1 : number_variables
                % The number of dimensions this variable occupies
                variable_name   = variable_names_cell{v};
                example_data    = Variable_Data(1).(variable_name);         % The first entry is used just as an example
                num_dim         = length(example_data);

                end_ind         = start_ind + num_dim - 1;
                
                % This variable's expected value, standard deviation and confidence interval
                expected_value              = variable_mu_list(start_ind : end_ind);
                standard_deviation          = variable_STD_list(start_ind : end_ind);
                confidence_interval_matrix  = variable_bounds_matrix(:, start_ind : end_ind);
                
                Variable_Uncertainty.(variable_name) = struct('mu', expected_value, 'sigma', standard_deviation,'confidence_interval', confidence_interval_matrix);
                
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

