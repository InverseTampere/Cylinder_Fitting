% This is an updated version of the script Pasi Raumonen wrote for TreeQSM

% --------------------------------------------------------------------------
% Input 
%   alpha       Rotation around x-axis 
%   beta        Rotation around y-axis

% Output 
% R             Full rotation matrix
% R_alpha       Rotation of alpha [1 0 0; 0 c1 -s1; 0 s1 c1]
% R_beta        Rotation of beta [c2 0 s2; 0 1 0; -s2 0 c2]
% dR_alpha      Derivative of R w.r.t. alpha
% dR_beta       Derivative of R w.r.t. beta
% ---------------------------------------------------------------------

%%%%%% CHANGES %%%%%%
% Changed the inputs to specifically be alpha and beta
% Changed the outputs to include the alpha and beta rotation matrices
% Wrote out the sine and cosine arguments

function [R, R_alpha, R_beta, dR_alpha, dR_beta] = form_rotation_matrices_2(alpha, beta)
  
    %% Rotation matrices %%
        % Alpha rotation matrix
        R_alpha = [1, 0, 0; 0 cos(alpha), -sin(alpha); 0, sin(alpha), cos(alpha)];

        % Beta rotation matrix
        R_beta  = [cos(beta), 0, sin(beta); 0, 1, 0; -sin(beta), 0, cos(beta)];

        % Full rotation matrix
        R       = R_beta * R_alpha;
    
    %% Their derivatives %%
        if nargout > 3
            dR_alpha    = [0, 0, 0; 0, -R_alpha(3, 2), -R_alpha(2, 2); 0, R_alpha(2, 2), -R_alpha(3, 2)];
        end
        
        if nargout > 4
            dR_beta     = [-R_beta(1, 3), 0, R_beta(1, 1); 0, 0, 0; -R_beta(1, 1), 0, -R_beta(1, 3)];
        end
end