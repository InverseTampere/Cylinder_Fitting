% Iteratively sampling the function to find its minima is a more robust, but potentially slower and more memory-intensive, alternative to golden-section search
% The function handle must have only one input, which is the variable over which is searched
% Note that the given convergence threshold is absolute and in terms of the function value
% The local minima are given in ascending order

function [global_min_var_value, global_min_fun_value, local_min_var_value_list, local_min_fun_value_list, number_local_minima] = Sampling_Function_Minimiser(Function_handle, search_lower_bound, search_upper_bound, convergence_threshold, max_iterations, Print, Plot)
    
    %% Manual inputs %%
        number_samples      = 1e2;      % [-] The number of samples during each iteration
        numerical_margin    = 1e-6;     % [-] To check whether the search space bounds were hit

    %% Detected local minima %%
        % Recursive function to find the local minima
        iter = 0;           % Initialising the iteration counter
        [local_min_var_value_list, local_min_fun_value_list, number_local_minima] = Function_Minimiser(Function_handle, search_lower_bound, search_upper_bound, iter, number_samples, convergence_threshold, max_iterations, Print);
        
        % They are sorted by magnitude
        [local_min_fun_value_list, order]   = sort(local_min_fun_value_list, 'ascend');
        local_min_var_value_list            = local_min_var_value_list(order);

        % The global minimum is now the first entry
        global_min_fun_value = local_min_fun_value_list(1);
        global_min_var_value = local_min_var_value_list(1);

        if Print == true
            fprintf('The global function minimum is %.3g at variable value %.3g \n', global_min_fun_value, global_min_var_value);
        end

        % A warning is displayed if Print is true and any of the variable values are at the edge of the search space
        if Print == true
            % Upper bound
            UB_diff_list = abs(local_min_var_value_list - search_upper_bound);
            
            if min(UB_diff_list) < numerical_margin                             % A margin is used
                warning('The upper bound of the search space was hit.')
            end
    
            % Lower bound
            LB_diff_list = abs(local_min_var_value_list - search_lower_bound);

            if min(LB_diff_list) < numerical_margin                             % A margin is used
                warning('The lower bound of the search space was hit.')
            end
        end

    %% Plot %%
        % Plot to see whether the minima were found correctly
        if Plot == true
            % The function is sampled within the search interval
            search_interval_list    = linspace(search_lower_bound, search_upper_bound, 2*number_samples);       % Twice the number of samples is used now
            function_value_list     = arrayfun(Function_handle, search_interval_list);                

            % The y-axis of the plot is limited by the function values
            plot_upper_limit = max(function_value_list) + 0.1*(max(function_value_list) - min(function_value_list));
            plot_lower_limit = min(function_value_list) - 0.1*(max(function_value_list) - min(function_value_list));

            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   

            hold on
            grid on

            % The function values
            plot(search_interval_list, function_value_list, 'color', 'k', 'LineWidth', 2, 'DisplayName', 'Function values');

            % The found minima
            for m = 1 : number_local_minima
                variable_value  = local_min_var_value_list(m);
                pl_minimum      = plot(variable_value*[1,1], [plot_lower_limit, plot_upper_limit], 'color', 'm', 'LineWidth', 2, 'DisplayName', 'Found minimum');

                if m > 1
                    pl_minimum.HandleVisibility = 'off';
                end
            end

            % Axes
            xlabel('t');
            ylabel('f(t)');

            xlim([search_lower_bound, search_upper_bound]);
            ylim([plot_lower_limit, plot_upper_limit]);

            % Legend
            legend('show', 'location', 'northoutside');
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            % Paused for analysis
            disp('The local minima have been found. The figure will close and the script will end when a key is pressed.');
            pause();
            close(1);
        end
            
    %% Local iterative minimiser function %%
        % A local function is created s.t. it can be treated recursively
        function [local_minima_var_value_list, local_minima_fun_value_list, number_local_minima] = Function_Minimiser(Function_handle, search_lower_bound, search_upper_bound, iter, number_samples, convergence_threshold, max_iterations, Print)
            % Iteration counter is updated
            iter = iter + 1;

            % Function values between the search bounds
            t_list = linspace(search_lower_bound, search_upper_bound, number_samples);
            F_list = arrayfun(Function_handle, t_list);

            % The difference between each step
            dF_list = diff(F_list);

            % Local minima have a negative and positive derivative before and after respectively
            local_minima_ind    = find(dF_list(1 : number_samples - 2) < 0 & dF_list(2 : number_samples - 1) > 0);
            number_local_minima = length(local_minima_ind);

            if number_local_minima == 0                     % At the edge of the domain the derivative check may not hold
                [~, local_minima_ind]   = min(F_list);
                number_local_minima     =  1;
            end
            
            % Each local minimum is iterated upon if needed
            local_minima_var_value_cell = cell(1, number_local_minima);
            local_minima_fun_value_cell = cell(1, number_local_minima);

            for i = 1 : number_local_minima
                % Indices of the local minimum and its neighbours
                local_minimum_ind   = local_minima_ind(i);
                prev_ind            = max(1, local_minimum_ind - 1);
                next_ind            = min(number_samples, local_minimum_ind + 1);

                % Function values at the previous and next indices are checked for convergence
                F_prev      = F_list(prev_ind);
                F_next      = F_list(next_ind);
                fun_diff    = abs(F_prev - F_next);             % Note that this is a conservative convergence criterion

                if fun_diff < convergence_threshold             % The difference is small enough for convergence 
                    % The local minimum is appended
                    local_minima_var_value_cell{i} = t_list(local_minimum_ind);
                    local_minima_fun_value_cell{i} = F_list(local_minimum_ind);

                elseif iter == max_iterations                               % The maximum number of iterations was met
                    % The local minimum is appended insofar as it was found
                    local_minima_var_value_cell{i} = t_list(local_minimum_ind);
                    local_minima_fun_value_cell{i} = F_list(local_minimum_ind);

                    % A warning is displayed if Print is true
                    if Print == true
                        warning('The local minimum could not be found in %i iterations. The difference is still %.3g \n', fun_diff);
                    end
                else                                                        % Another iteration is required
                    % The search bounds are updated
                    search_lower_bound_i = t_list(prev_ind);
                    search_upper_bound_i = t_list(next_ind);
                    
                    % The local minima are appended to the cell array
                    [local_minima_var_value_list_i, local_minima_fun_value_list_i, ~] = Function_Minimiser(Function_handle, search_lower_bound_i, search_upper_bound_i, iter, number_samples, convergence_threshold, max_iterations, Print);
                    local_minima_var_value_cell{i} = local_minima_var_value_list_i;
                    local_minima_fun_value_cell{i} = local_minima_fun_value_list_i;
                end
            end

            % All found local minima
            local_minima_var_value_list = horzcat(local_minima_var_value_cell{:});
            local_minima_fun_value_list = horzcat(local_minima_fun_value_cell{:});
            number_local_minima         = length(local_minima_fun_value_list);
        end
end
