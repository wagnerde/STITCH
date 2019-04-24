function  [unique_object_counts, unique_objects] = count_unique(list, plot_hist)
%% Usage: [unique_object_counts, unique_objects] = count_unique(list, plot_hist)
%
% Based on the built-in Matlab function: unqiue.  Returns unique elements
% from a list along with their respective counts.
%

%% SETTINGS:
if ~exist('plot_hist', 'var')
    plot_hist = false;
end

%% CODE:
[unique_objects, ~, ic] = unique(list);
unique_object_ind = 1:length(unique_objects);
unique_object_counts = histc(ic,unique_object_ind)';

%% OPTIONAL HISTOGRAM PLOT
if plot_hist
    figure
    bar(1:length(unique_objects), unique_object_counts)
    ylabel('Abundance')
    xlabel('Unique Objects')
end