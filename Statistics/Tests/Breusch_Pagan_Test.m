% The Breusch-Pagan test (Breusch & Pagan, 1979) can be used to find out whether the given data is heteroscedastic or not
% Optionally, it may be studentised if the error term is non-Gaussian (Koenker, 1981) 

% Regression may be performed beforehand or by the script, requiring two different input structures

% Regression performed here:
    % Regression_type       ('Polynomial' or 'Circle')
    % polynomial_order      (Only required for a polynomial) 
    % predictor_list        (x for polynomial, theta for circle)
    % response_list         (y for polynomial, r for circle)

% Beforehand:
    % residual_list         (List of residuals for each predictor value)
    % predictor_list        (i.e. x for a polynomial)
    % Degrees_of_Freedom    (i.e. polynomial order)

    function [BP_statistic, P_value] = Breusch_Pagan_Test(Input_Structure, Studentise)
    
    %% Regression %%
        % Regression was already performed if the residuals were already given
        try
            % Data from the input structure
            residual_list       = Input_Structure.residual_list;
            predictor_list      = Input_Structure.predictor_list;
            Degrees_of_Freedom  = Input_Structure.Degrees_of_Freedom;

            % Ensured to be vertical
            n               = length(predictor_list);
            predictor_list  = reshape(predictor_list, [n, 1]);
            residual_list   = reshape(residual_list, [n, 1]);
    
        catch
            % Data from the input structure
            Regression_type     = Input_Structure.Regression_type;
            predictor_list      = Input_Structure.predictor_list;
            response_list       = Input_Structure.response_list;

            % Ensured to be vertical
            n               = length(predictor_list);
            predictor_list  = reshape(predictor_list, [n, 1]);
            response_list   = reshape(response_list, [n, 1]);
    
            if strcmp(Regression_type, 'Polynomial')
                % Polynomial order
                polynomial_order    = Input_Structure.polynomial_order;
                Degrees_of_Freedom  = polynomial_order;

                % Matrix of the predictor 
                predictor_matrix    = predictor_list.^(1 : polynomial_order);
    
                % Fitting the model
                Regression_Model    = fitlm(predictor_matrix, response_list);
    
                % Residuals
                residual_list       = Regression_Model.Residuals.Raw;
    
            elseif strcmp(Regression_type, 'Circle')
                % The degrees of freedom is one
                Degrees_of_Freedom = 1;

                % Least-squares circle fitting
                point_matrix = response_list .* [cos(predictor_list), sin(predictor_list)];
    
                A = [point_matrix, ones(n, 1)];
                B = -sum(point_matrix.^2, 2);
                
                circle_coefficients = A \ B;
                circle_centre       = -1/2 * circle_coefficients(1:2)'; 
                circle_radius       = sqrt(sum(circle_centre.^2) - circle_coefficients(3));
    
                % Residuals
                point_matrix_c  = point_matrix - circle_centre;
    
                residual_list   = sum(point_matrix_c.^2, 2) - circle_radius^2;
    
            else
                error('The regression type must be Polynomial or Circle.');
            end
        end

    %% Breusch-Pagan test %%
        % Fitting the auxiliary regression model to the residuals
        squared_residuals   = residual_list.^2;
        Auxiliary_Model     = fitlm(predictor_list, squared_residuals);

        % R-squared value
        R_squared_value = Auxiliary_Model.Rsquared.Ordinary;

        % Studentised Breusch-Pagan test statistic
        n               = length(residual_list);
        BP_statistic    = R_squared_value * n;

        % Adjusted if a non-studentised statistic is desired
        if ~Studentise
            adjustment_factor   = (n - 1)/n * var(residual_list.^2) / (2*((n - 1)/n * var(residual_list)).^2);
            BP_statistic        = adjustment_factor * BP_statistic;
        end

        % P-value of the Breusch-Pagan test
        P_value = 1 - chi2cdf(abs(BP_statistic), Degrees_of_Freedom);
end








