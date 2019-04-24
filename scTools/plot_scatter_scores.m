function plot_scatter_scores(xy, scores)
% USAGE: plot_scatter_scores(xy, scores)
%
% Generates a simple 2d scatterplot for overlaying continuous scoring data
%
% INPUTS:
% xy            2d array of xy coordinates
% scores        Numerical scores for each row in xy
%

%% CODE:

% make scatter plot
figure('Position', [500 500 500 500])
scatter(xy(:,1), xy(:,2), 15, [0.7 0.7 0.7], 'fill','linewidth',0.0001,'MarkerFaceAlpha',0.8);
hold on
scatter(xy(:,1), xy(:,2), 25, scores, 'fill','linewidth',0.0001,'MarkerFaceAlpha',0.8);

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
colormap(jet)



