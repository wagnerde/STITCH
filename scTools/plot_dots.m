function plot_dots(X, gene_names, gene_names_all, cell_ind_groups, group_names)
% Usage: plot_dots(X, gene_names, gene_names_all, cell_ind_groups, group_names)
%
% Generates a "Dot Plot" style heatmap summarizing gene expression patterns
% over a series of cell groups. In the heatmap, dot diameter indicates
% the fraction of cells in each group that are positive for each gene 
% (i.e. contain non-zero counts).  Dot color indicates mean expression
% value for each gene.
% 
% INPUTS:
% 'X' 
%       A total-count-normalized transcripts x cells counts matrix 
%       (rows=transcripts; columns=cells).
%
% 'gene_names'
%       Cell array of strings. Names of genes to be displayed as rows in 
%       the heatmap.
%
% 'gene_names_all'
%       Cell array of strings. Names for all genes (corresponding to rows in X).
%
% 'cell_ind_groups'
%       Cell array containing a vector of cell indices. Each vector 
%       corresponds to a cell group, and its values indicate column indices
%       in X for each cell contained in that group.
% 
% 'group_names'
%       Cell array of strings.  Labels for each cell group organized as
%       columns in the heatmap.
%

%% CODE:

% convert gene names to row indices
gene_names = fliplr(gene_names);
for j = 1:length(gene_names)
    [~, gene_ind(j)] = ismember(gene_names{j}, gene_names_all);
end

% get mean expression and fraction expressing cells for each gene / group 
nColumns = length(cell_ind_groups);
for k = 1:nColumns  
    % subset by sample
    X_sub = X(gene_ind, cell_ind_groups{k});
    % calculate mean expression and fPositive cells for each gene
    Mean{k} = mean(X_sub,2);
    fPos{k} = sum(X_sub~=0,2)/size(X_sub,2);
end

% format and scale data
grid_fPos = cell2mat(fPos);
% scale gene expression by zscore
grid_Z = zscore(cell2mat(Mean)); 
% scale gene expression by range
grid_Range = (cell2mat(Mean) - min(cell2mat(Mean),[],2)) ./ (max(cell2mat(Mean),[],2) - min(cell2mat(Mean),[],2));

% make grid
nRows = length(gene_names);
[grid_X, grid_Y] = meshgrid(1:nColumns,1:nRows);

% make figure
figure('Position', [200 200 300 600])
set(gcf, 'Color', 'white')
ax1 = axes('Position', [0.1 0.1 0.3 0.8]);
scatter(grid_X(:), grid_Y(:), 'Filled', 'CData', grid_Range(:), 'SizeData', 5+95*grid_fPos(:))
set(ax1,'YAxisLocation', 'right', 'XAxisLocation', 'top');
set(ax1, 'XTick', (1:nColumns), 'XTickLabel', group_names, 'YTick', (1:nRows), 'YTickLabel', gene_names, 'TickLength',[0 0]) 
xtickangle(90)
box on
xlim([0 nColumns+1])
ylim([0 nRows+1])
clim = get(gca, 'Clim');
col = jet(1000);
colormap(col(1:920,:)) % get rid of the brown

% colorbar legend
ax2 = axes('Position', [0.4 0.1 0.4 0.18]);
c = colorbar('Location', 'eastoutside');
caxis(clim)
c.Label.String = 'Relative Expression';
axis off

% fPos legend
ax3 = axes('Position', [0.66 0.35 0.06 0.1]);
c = colorbar('Location', 'eastoutside');
scatter([1 1 1],[1 2 3],'Filled','CData',[0 0 0], 'SizeData', 5+95*[0.1 0.5 1]); 
set(ax3, 'YTick', [1 2 3], 'YTickLabel', {'10%','50%','100%'}', 'XTickLabel', [],'TickLength',[0 0]) 
set(ax3,'YAxisLocation', 'right')
text(1.2, 4.5, {'Fraction','Positive Cells'}, 'FontWeight', 'bold', 'HorizontalAlignment', 'center'); text(1.15, 3, '100%'); text(1.15, 2, '50%'); text(1.15, 1, '10%'); 
xlim([0.9 1.1])
ylim([0 4])
caxis(clim)
axis off

