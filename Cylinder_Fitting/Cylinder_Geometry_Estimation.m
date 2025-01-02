% The cylinder geometry is determined in two stages using likelihood maximisation (and distance minimisation)
% First the cross-section and direction are determined, then the length

function [Optimal_Geometry, LS_Geometry, Optimiser_Diagnostics, Point_Cloud_Distributions] = Cylinder_Geometry_Estimation(Point_Cloud_Coord, Scanning_Parameters, Scanner_Parameters, Statistical_Values, Fitting_Parameters, Output_Decisions)

    %% Cross-section and direction estimation %%
        % This geometry is estimated through likelihood and distance optimisation
        % This geometry is also used to estimate the uncertainty
        [Optimal_Geometry, LS_Geometry, ~, Point_Cloud_Distributions, Optimiser_Diagnostics] = Cylinder_Cross_Section_and_Direction_Estimation_fmincon(Point_Cloud_Coord, Scanning_Parameters, Scanner_Parameters, Statistical_Values, Fitting_Parameters, Output_Decisions);
        
        % The length is estimated by maximising the likelihood of the given point cloud
        Cylinder_centre     = Optimal_Geometry.Cylinder_centre;
        Cylinder_direction  = Optimal_Geometry.Cylinder_direction;
        [Cylinder_length, Cylinder_top_loc, Cylinder_bot_loc] = Cylinder_Length_Estimation(Cylinder_centre, Cylinder_direction, Fitting_Parameters, Point_Cloud_Distributions, Output_Decisions);

        % The length parameters are added to the optimal geometry structure
        Optimal_Geometry.Cylinder_length    = Cylinder_length;
        Optimal_Geometry.Top_loc            = Cylinder_top_loc;
        Optimal_Geometry.Bottom_loc         = Cylinder_bot_loc;
        
end