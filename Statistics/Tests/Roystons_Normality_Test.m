% This script uses Royston's test to determine whether the data is normally distributed or not at confidence level alpha
% Normal distribution is indicated by the Normal_Flag being true

% Note: The expected matrix is expected to be of shape [N, M], where M is the number of variables
% Note: N needs to be between 3 and 2000

function [Roystons_Statistic, P_value, Normal_Flag] = Roystons_Normality_Test(Data_Matrix, alpha, Print)

    %% Shapiro-Francia / Shapiro-Wilk's test %%
        % The amount of given data
        [N, M] = size(Data_Matrix);

        if N < 3 || N > 2000
            error('Error: The Royston''s test requires the number of observations to be between 3 and 2000.');
        end
        
        % The test statistic values
        W_list = zeros(1, M);
        
        for j = 1 : M
            % This column's sorted (centered) data
            x_list      = Data_Matrix(:, j);
            x_list      = sort(x_list);
            x_list_c    = x_list - mean(x_list);
            
            m_list = norminv(((1 : N) - 3/8)/(N + 1/4))';
            c_list = m_list / sqrt(m_list' * m_list);
            
            % The kurtosis dictates which test is performed
            K = kurtosis(x_list);
            
            if K > 3
                %--% The Shapiro-Francia test %--%                
                W = (c_list' * x_list)^2 / (x_list_c' * x_list_c);
                
            else
                %--% The Shapiro-Wilk's test %--%
                % Polynomial coefficients
                P_1 = [-2.706056, 4.434685, -2.071190, -0.147981, 0.221157, c_list(N)];
                P_2 = [-3.582633, 5.682633, -1.752461, -0.293762, 0.042981, c_list(N - 1)];            
                
                w_list = zeros(N, 1);
                
                if N == 3
                    w_list(1) = 1/sqrt(2);
                    w_list(N) = -w_list(1);
                else
                    w_list(N) = polyval(P_1, 1/sqrt(N));
                    w_list(1) = -w_list(N);
                end
                
                if N < 6
                    t   = 2;
                    phi = (m_list' * m_list - 2*m_list(N)^2) / (1 - 2*w_list(N)^2);
                else
                    w_list(N - 1)   = polyval(P_2, 1/sqrt(N));
                    w_list(2)       = -w_list(N - 1);
                    
                    t   = 3;
                    phi = (m_list' * m_list - 2*(m_list(N)^2 + m_list(N - 1)^2)) / (1 - 2*(w_list(N)^2 + w_list(N - 1)^2));                    
                end
            
                w_list(t : N - t + 1) = m_list(t : N - t + 1) / sqrt(phi);
                
                % The test statistic
                W = (w_list' * x_list)^2 / (x_list_c' * x_list_c);
            end
            
            % The test statistic is appended            
            W_list(j) = W;
        end

    %% Royston's test %%
        if N < 12
            g = -2.273 + 0.459*N;
            m = 0.5440 - 0.39978*N + 0.025054*N^2 - 0.0006714*N^3;
            s = exp(1.3822 - 0.77857*N + 0.062767*N^2 - 0.0020322*N^3); 
            
            Z_list = 1/s * (-log(g - log(1 - W_list)) - m);
        else
            m = -1.5861 - 0.31082*log(N) - 0.083751*log(N)^2 + 0.0038915*log(N)^3;
            s = exp(-0.4803 -0.082676*log(N) + 0.0030302*log(N)^2);
            
            Z_list = 1/s * (log(1 - W_list) - m);
        end
    
        % Correlation of the data matrix
        C = corrcoef(Data_Matrix);
        
        u = 0.715;
        v = 0.21364 + 0.015124*(log(N))^2 - 0.0018034*(log(N))^3;
        l = 5;

        % Transformed correlation matrix
        C_T = (C.^l) .* (1 - (u*(1 - C).^u) / v);
        
        % The average correlation
        C_total = sum(sum(C_T)) - M;
        c_avg   = C_total / (M^2 - M);
        
        % The degrees of freedom
        if M > 1
            DoF = M/(1 + (M - 1) * c_avg);
        else
            DoF = 1;
        end
        
        % Royston's statistic value
        R_list              = norminv(normcdf(-Z_list) / 2).^2;
        
        Roystons_Statistic  = DoF / M * sum(R_list);
        
        % The statistical significance
        P_value = 1 - chi2cdf(Roystons_Statistic, DoF);
        
        if P_value >= alpha
            Normal_Flag = true;
        else
            Normal_Flag = false;
        end
        
        % Printed statements
        if Print == true
            fprintf('\n')
            disp('--------------------')
            if Normal_Flag == true
                fprintf('The given data matrix is normally distributed, with P = %.3g >= alpha = %.2g \n', P_value, alpha);
            else
                fprintf('The given data matrix is NOT normally distributed, with P = %.3g < alpha = %.2g \n', P_value, alpha);
            end
            
            fprintf('The data has %g equivalent degrees of freedom \n', DoF);
        end
end