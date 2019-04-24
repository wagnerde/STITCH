function plot_scatter_discrete(xy, groupIDs)
% USAGE: plot_scatter_discrete(xy, groupIDs)
%
% Generates a simple 2d scatterplot for visualizing discrete cell groups 
% (e.g. clusters).
%
% INPUTS:
% xy            2d array of xy coordinates
% groupIDs      Group assignments for each cell / each row in xy
%               Can be numerical, or strings
%

%% CODE:

% if group names are not specified, just use numbers instead
if iscell(groupIDs)
    group_ind = grp2idx(groupIDs);
    group_names = groupIDs;
elseif isnumeric(groupIDs)
    group_ind = groupIDs;
    group_names = cellstr(num2str(groupIDs));
end

% use cbrewer to generate color heatmaps
% https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
addpath('scTools/cbrewer/') 

% set colormap
cols = cbrewer('qual','Set3', max(group_ind));

% if 'unassigned' (i.e. "-1") ID is present, set to black
if sum(ismember(group_ind, -1))>0
    cols = [cols; [0 0 0]];
    group_ind(group_ind == -1) = max(group_ind)+1;
end

% make scatter plot
%figure('Position', [500 500 500 500])
scatter(xy(:,1), xy(:,2), 25, group_ind, 'fill','linewidth',0.0001,'MarkerFaceAlpha',0.8);

% adjust axis limits to leave gap
gap = 0.1;
x_plot = xy(:,1); x_range = range(x_plot);
y_plot = xy(:,2); y_range = range(y_plot);
set(gca,'xlim', [min(x_plot)-gap*x_range, max(x_plot)+gap*x_range])
set(gca,'ylim', [min(y_plot)-gap*y_range, max(y_plot)+gap*y_range])    

% adjust plot appearance
axis square; axis off;
set(gca,'xticklabel',{''},'yticklabel',{''});
set(gcf,'color','w');
colormap(cols)

% add group labels
for k = 1:max(group_ind)
    text(mean(xy(group_ind==k,1)), mean(xy(group_ind==k,2)), group_names(find(group_ind==k,1)), 'fontweight', 'normal', 'fontsize', 12, 'Interpreter', 'none');
end


