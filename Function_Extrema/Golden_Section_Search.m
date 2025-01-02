% This script finds the minimum of the given function using golden-section search
% The function handle must have only one input, which is the variable over which is searched
% Note that the given convergence threshold is absolute and in terms of the function value

function [variable_value, function_minimum, final_search_bounds] = Golden_Section_Search(Function_handle, search_lower_bound, search_upper_bound, convergence_threshold, max_iterations, Print, Diagnostics)
    
    %% Manual inputs %%
        numerical_margin    = 1e-6;     % [-] Small value to determine whether the end result lies on either of the bounds

    %% Golden-section search %%
        % The golden ratio
        psi = (1 + sqrt(5))/2;

        % Iteration prep
        convergence = false;
        iter        = 0;
        F_old       = Inf;

        t1 = search_lower_bound;
        t4 = search_upper_bound;
        
        while convergence == false && iter < max_iterations
            iter = iter + 1;

            % The points at the golden ratios from t1 and t4
            t2 = t1 + (2 - psi) * (t4 - t1);
            t3 = t1 + (psi - 1) * (t4 - t1);

            % Function values at t1 to t4
            t_list  = [t1, t2, t3, t4];
            
            F_list = arrayfun(Function_handle, t_list);

            % The lowest two points
            [F_min_list, min_ind_list]  = sort(F_list, 'ascend');
            F_min                       = F_min_list(1);
            [min_ind_1, min_ind_2]      = deal(min_ind_list(1), min_ind_list(2));
            
            % Convergence check
            if abs(F_min - F_old) < convergence_threshold
                convergence = true;
                
            % Without convergence, the points are updated
            else
                F_old = F_min;
                
                t1 = t_list(min([min_ind_1, min_ind_2]));
                t4 = t_list(max([min_ind_1, min_ind_2]));
            end
        end
        
    %% Result %%
        % The minimum function, respective variable value and the final search interval
        variable_value      = t_list(min_ind_1);
        function_minimum    = Function_handle(variable_value);
        final_search_bounds = [t1, t4];
        
        if Print == true
            fprintf('After %g iterations the minimal function value is %.3g, with variable value %.3g \n', iter, function_minimum, variable_value);
        end
        
    %% Checks %%
        % An error is printed if the minimum value is NaN
        if isnan(F_min)
            error('The found function value is NaN')
        end
        
        % Warnings are only printed if Print is true
        if Print == true
            num_min_entries = length(find(F_list == F_min));

            if num_min_entries > 1
                warning('The minimum function value was found at multiple locations');
            end

            % If the maximum number of iterations was hit, a message is printed
            if iter == max_iterations
                warning('The maximum number of iterations, %g, was reached', max_iterations);
            end

            % If the final variable value is at the edge of the initial interval, a different warning is printed
            diff_LB = abs(variable_value - search_lower_bound);
            diff_UB = abs(variable_value - search_upper_bound);

            if diff_LB < numerical_margin || diff_UB < numerical_margin      % A slight margin is used in case of rounding errors
                warning('The minimum value was at the edge of the search interval');
            end
        end
        
        % Diagnostics plot to see whether the minimal function value was found correctly
        if Diagnostics == true
            % The function is sampled within the search interval
            number_samples          = 1e3;
            search_interval_list    = linspace(search_lower_bound, search_upper_bound, number_samples);
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
            
            % The final sample locations
            for i = 1 : 4
                t       = t_list(i);
                pl_t    = plot(t*[1, 1], [plot_lower_limit, plot_upper_limit], 'color', 'b', 'LineWidth', 2, 'DisplayName', 'Final sample locations');
                    
                if i > 1
                    pl_t.HandleVisibility = 'Off';
                end
            end
            
            % The found and true minima
            plot(variable_value*[1,1], [plot_lower_limit, plot_upper_limit], 'color', 'm', 'LineWidth', 2, 'DisplayName', 'Found minimum');
            
            [~, min_ind]    = min(function_value_list);
            true_minimum    = search_interval_list(min_ind); 
            plot(true_minimum*[1,1], [plot_lower_limit, plot_upper_limit], 'color', 'c', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'True minimum');
            
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
            disp('The figure will close and the script will continue when a key is pressed.');
            pause();
            close(1);
        end
        
end