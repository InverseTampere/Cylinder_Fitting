% This script determines the optimal number of bins using the Freedman-Diaconis rule

function [number_bins, bin_width] = Histogram_Bins(histogram_data)

    % Data properties
    IQR             = iqr(histogram_data);                              % Interquartile range
    data_amplitude  = max(histogram_data) - min(histogram_data);

    if IQR < 1e-12                                                      % If there is no variety in the data, it is placed in a single bin
        number_bins = 1;
        bin_width   = data_amplitude;
    else
        bin_width       = 2 * IQR * length(histogram_data)^(-1/3);      % Freedman-Diaconis rule
        number_bins     = round(data_amplitude / bin_width);
        number_bins     = max(1, number_bins);                          % In case the number of samples is 1
    end

    % The number of bins can become excessive if there are outliers in the data. More than 1000 bins are difficult to distinguish graphically, so the limit is set there
    number_bins = min(number_bins, 1e3);

end