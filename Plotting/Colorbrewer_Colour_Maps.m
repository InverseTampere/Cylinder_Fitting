% This script returns multiple unique diverging('Diverging' / 'div') or sequential('Sequential' / 'seq') colour maps with the given number of colours in the given category
% The colour map categories are the same as in cbrewer: 'div', 'seq'

function colour_map_cell = Colorbrewer_Colour_Maps(colour_map_category, number_colour_maps, number_colours)

    %% Colour map names %%
        % Diverging. Note that only twin colour maps are used
        div_cmap_name_cell  = {'RdGy', 'RdBu', 'PuOr', 'PRGn', 'PiYG', 'BrBG'};

        % Sequential. All colour maps are used, but the order is changed
        seq_cmap_name_cell  = {'Reds', 'Purples', 'Greens', 'Blues', 'Oranges', 'Greys', 'YlOrRd', 'YlGnBu', 'YlGn', 'RdPu', 'PuRd', 'PuBuGn', 'OrRd', 'GnBu', 'BuPu', 'BuGn', 'PuBu', 'YlOrBr'};            

    %% Colour maps %%
        % Selecting the right data
        if strcmp(colour_map_category, 'Diverging') || strcmp(colour_map_category, 'div')
            cbrewer_category    = 'div';
            cmap_name_cell      = div_cmap_name_cell;
        elseif strcmp(colour_map_category, 'Sequential') || strcmp(colour_map_category, 'seq')
            cbrewer_category    = 'seq';
            cmap_name_cell      = seq_cmap_name_cell;
        else
            error('The colour map category must be Diverging or Sequential');
        end

        % Repeating the colour maps if needed
        number_unique_cmaps = length(cmap_name_cell);
        
        if number_colour_maps > number_unique_cmaps
            warning('There aren''t enough unique colour maps, so they need to be repeated.');

            number_repeats          = ceil(number_colour_maps / number_unique_cmaps);
            cmap_name_cell_total    = cell(1, number_repeats);
            cmap_name_cell_total(:) = {cmap_name_cell};
            cmap_name_cell          = horzcat(cmap_name_cell_total{:});
        end

        % Generating the colour maps
        colour_map_cell = cell(1, number_colour_maps);

        for c = 1 : number_colour_maps
            % This colour map is generated
            cmap_name       = cmap_name_cell{c};
            cmap_matrix     = cbrewer(cbrewer_category, cmap_name, max(3, number_colours));         % Note that a minimum of 3 is required
            cmap_matrix     = cmap_matrix(1 : number_colours, :);                                   % In case there were fewer than 3, the excess is removed

            % Sometimes interpolation can lead to values outside of the normalised bounds
            cmap_matrix     = max(cmap_matrix, 0);
            cmap_matrix     = min(cmap_matrix, 1);
            
            colour_map_cell{c} = cmap_matrix;
        end
end