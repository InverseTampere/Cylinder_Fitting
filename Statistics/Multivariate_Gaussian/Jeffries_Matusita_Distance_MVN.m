% The Jeffries-Matusita (and thus also Bhattacharyya) distance is computed between two multivariate normal distributions
% Note that mu is expected to be a 1 x m vector

function [JM_distance, B_distance, B_coeff] = Jeffries_Matusita_Distance_MVN(mu_a, Sigma_a, mu_b, Sigma_b)

    %% Bhattacharyya distance %%
        % The average covariance matrix
        Sigma = (Sigma_a + Sigma_b) / 2;
    
        % The Bhattacharyya distance and coefficient
        B_distance  = 1/8*((mu_a - mu_b) / Sigma) * (mu_a - mu_b)' + 1/2*log(det(Sigma)/sqrt(det(Sigma_a)*det(Sigma_b)));            
        B_coeff     = exp(-B_distance);

    %% Jeffries-Matusita distance %%
        % It is now bounded between 0 and 2
        JM_distance = 2*(1 - B_coeff);
        
end
