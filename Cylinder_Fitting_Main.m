% The cylinder is fitted in two stagess. First the cross-section and direction are estimated, then the length

beep off
clear variables
close all
clc 

%% Inputs %%
    %--% Execution %--%
    rng(1);                                             % [-] Integer below 100. RNG can be set for testing purposes. Only affects the initial point cloud uncertainty

    Parallel_Loop           =   true;                  % [true, false] Determines whether or not the Monte Carlo loop is ran in parallel
    idle_timeout            =   030;                    % [min] Time until the parallel pool shuts down if idle        

    %--% Scanning parameters %--%
    % Scanning
    Initial_Uncertainty     =   false;                  % [true, false] Whether or not the input point cloud to the Monte Carlo loop is uncertain or not
    Scanner_loc_cell        =   {[0, 50, 0]};           % [m] Location of the scanners {[x1, y1, z1], [x2, y2, z2], ...} 

    % Cylinder geometry
    Cylinder_centre         =   [0.0, 0.0, 0.0];        % [m] Location of the centre of the cylinder (x, y, z)
    Cylinder_direction      =   [0.0, 0.0, 1.0];        % [-] The direction vector of the cylinder         
    Cylinder_length         =   0.25;                   % [m] Length of the cylinder (symmetrical around the centre)
    Cylinder_radius         =   0.05;                   % [m] Cylinder radius

    %--% Fitting parameters %--%    
    % Least-squares (initial) cylinder fitting parameters
    LS_Filtering                =   false;              % [true, false] Whether or not points outside the confidence interval are removed when LS fitting is performed
    RANSAC                      =   false;              % [true, false] Whether or not RANSAC is used for the initial cylinder fit

    % Infinite cylinder fitting parameters
    distance_moment             =   2;                  % [1, 2] The moment of the expected distance
    Point_Weighting             =   'None';             % [None, Info, Reliability, Value] The points are weighted over the whole loop by the set criterium
    bounds_margin               =   1.00;               % [-] Factor by which the geometry parameters can deviate from the initial estimates [-]
    Distance_Computation        =   'Line_Approx';      % [Numerical, Line_Approx] The expected Mahalanobis distance to the circular cross-section can be computed numerically 
                                                        %                          or analytically by approximating the circle as a tangent-line (Line_Approx) which is much faster.
    max_UU_iterations           =   1e1;                % [-] Maximum number of times the geometry and uncertainty may be updated each time the cylinder geometry is estimated
    UU_convergence_threshold    =   1e-3;               % [-] Convergence threshold in terms of the normalised geometry vector when updating the geometry and uncertainty

    % Length estimation
    alignment_threshold         =   45;                 % [deg] Angle by which the normal vector of a point's triangle must be aligned with the cylinder axis for it to be an edge point candidate

    %--% Laser scanner parameters %--%
    laser_resolution        =   2*0.1534;               % [mrad] Angular resolution between pulses
    angular_accuracy        =   0.092;                  % [mrad] Accuracy of the angular encoder 
    beam_divergence         =   1*0.60;                 % [mrad] Divergence angle that covers 2 sigma of the beam power (half-angle)
    beam_exit_diameter      =   4.24;                   % [mm] Diameter of +/-2 sigma of the beam at the scanner exit
    range_bias              =   00.0;                   % [mm] The range bias
    sigma_range_device      =   01.0;                   % [mm] The range uncertainty of the instrument itself
    max_incidence_angle     =   80;                     % [deg] Max incidence angle can be limited to prevent points that are near-oblique, 0 to 90

    %--% Monte Carlo simulation %--%
    max_MC_length           =   050;                    % The maximum number of iterations of the Monte Carlo loop [-]
    number_GMM_fitment_iter =   1e2;                    % The maximum number of iterations for a Gaussian mixture model to be fitted [-]

    %--% Statistical %--%
    Confidence_interval     =   090;                    % [%] The percentage of the confidence interval. Note that this value is also used by RANSAC as the inlier threshold fraction
    Bhatt_coeff_threshold   =   Inf;                   % [-] Required Bhattacharyya coefficient for convergence (1 means complete probability overlap)

    %--% Outputs %--%
    Print                   =   true;                  % [true, false] Shows printed statements regarding intermediate and final results
    Plot                    =   true;                  % [true, false] Shows plots of final results

%% Input structures %%
    %--% Conversion to SI units %--%
    % Cylinder geometry
    Cylinder_direction  = Cylinder_direction / norm(Cylinder_direction);    % unit length

    % Fitting parameters
    alignment_threshold = deg2rad(alignment_threshold);                     % deg to rad

    % Scanner parameters
    laser_resolution    = laser_resolution * 1e-3;                          % mrad to rad
    angular_accuracy    = angular_accuracy * 1e-3;                          % mrad to rad
    beam_divergence     = beam_divergence * 1e-3;                           % mrad to rad
    beam_exit_diameter  = beam_exit_diameter * 1e-3;                        % mm to m
    range_bias          = range_bias * 1e-3;                                % mm to m
    sigma_range_device  = sigma_range_device * 1e-3;                        % mm to m
    max_incidence_angle = deg2rad(max_incidence_angle);                     % deg to rad

    %--% Creation of the structures %--%
    Parallel_Pool           = struct('Parallel_Loop', Parallel_Loop, 'idle_timeout', idle_timeout);
    True_Cylinder_Geometry  = struct('Cylinder_centre', Cylinder_centre, 'Cylinder_direction', Cylinder_direction, 'Cylinder_length', Cylinder_length, 'Cylinder_height_top', Cylinder_length/2, 'Cylinder_height_bot', -Cylinder_length/2, 'Cylinder_radius', Cylinder_radius);
    Scanning_Parameters     = struct('Initial_Uncertainty', Initial_Uncertainty, 'Scanner_loc_cell', {Scanner_loc_cell}, 'number_scanners', length(Scanner_loc_cell));
    Scanner_Parameters      = struct('laser_resolution', laser_resolution, 'beam_divergence', beam_divergence, 'beam_exit_diameter', beam_exit_diameter, 'range_bias', range_bias, 'sigma_range_device', sigma_range_device, 'angular_accuracy', angular_accuracy, 'max_incidence_angle', max_incidence_angle);
    Fitting_Parameters      = struct('LS_Filtering', LS_Filtering, 'RANSAC', RANSAC, 'distance_moment', distance_moment, 'bounds_margin', bounds_margin, 'Distance_Computation', Distance_Computation, 'Point_Weighting', Point_Weighting, 'alignment_threshold', alignment_threshold, 'max_UU_iterations', max_UU_iterations, 'UU_convergence_threshold', UU_convergence_threshold);
    Monte_Carlo_Inputs      = struct('max_MC_length', max_MC_length, 'number_GMM_fitment_iter', number_GMM_fitment_iter);
    Statistical_Values      = struct('Confidence_interval', Confidence_interval, 'Bhatt_coeff_threshold', Bhatt_coeff_threshold);
    Output_Decisions        = struct('Print', Print, 'Plot', Plot);

%% Point cloud %%
    tic;
    [Point_Cloud_Coord, Point_Cloud_Distributions] = Synthetic_Cylinder_Point_Generation(True_Cylinder_Geometry, Scanning_Parameters, Scanner_Parameters, Statistical_Values, Output_Decisions);
    t_Point_Cloud = toc;

    fprintf('   Generating the point cloud took %.5g seconds \n', t_Point_Cloud);

%% Optimal cylinder fit and uncertainty %%
    tic;
    [Optim_Cylinder_Geometry, Gaussian_Mixture_Models, Monte_Carlo_Loop, Fitting_Parameters] = Value_Weighted_Monte_Carlo_Loop_Cylinder(Parallel_Pool, Scanner_Parameters, Scanning_Parameters, Monte_Carlo_Inputs, Statistical_Values, Fitting_Parameters, Point_Cloud_Coord, Output_Decisions);
    t_MC = toc;

    fprintf('   Monte Carlo cylinder fitting took %.5g seconds \n', t_MC);

%% Result evaluation %% 
    tic;
    Cylinder_Result_Evaluation(True_Cylinder_Geometry, Point_Cloud_Coord, Optim_Cylinder_Geometry, Gaussian_Mixture_Models, Monte_Carlo_Loop, Scanning_Parameters, Scanner_Parameters, Fitting_Parameters, Statistical_Values, Output_Decisions);
    t_Results = toc;

    fprintf('   Result evaluation took %.5g seconds \n', t_Results);
