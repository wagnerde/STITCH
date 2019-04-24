function DynGenesTable = get_dynamic_genes(X, traj_pos, gene_names, sliding_window, counts_thresh, min_cells, fdr_alpha, gene_names_show)
% Usage: DynGenesTable = get_dynamic_genes(X, traj_pos, gene_names, sliding_window, counts_thresh, min_cells, fdr_alpha, gene_names_show)
%
% This function identifies genes expressed dynamically along a trajectory
% by implementing a method similar to that of Macosko et al (Cell 2015).  
% For each gene, a sliding window is first scanned along the trajectory to
% identify two windows with maximum and minimum average expression,
% respectively.  A t-test (fdr-corrected) is then used to identify
% significant expression level differences amongst individual cells from
% the two windows. Genes displaying significant differential expression are
% then ordered by peak expression along the trajectory.  
%
% INPUTS:
% 'X'
%       A total-count-normalized transcripts x cells counts matrix 
%       (rows=transcripts; columns=cells).
%
% 'traj_pos'
%       A numerical vector of trajectory positions for each cell in X 
%       (corresponding to columns in X).  NaN values will be ignored.  
% 
% 'gene_names'
%       Cell array of strings.  Names for all genes (corresponding to rows
%       in X).
%
% 'sliding_window'
%       Size of the sliding window, expressed as # of cells (e.g. '100').
%
% 'counts_thresh' 
%       Filter genes based on minimum number of counts (e.g. '3').
%
% 'min_cells' 
%       Minimum number of cells that must exceed counts_thresh (e.g. '1').
%
% 'fdr_alpha' 
%       Significance threshold for fdr-corrected p-value (e.g. '0.05').
%
% 'gene_names_show'
%       Optional. Names of genes to display on the final heatmap.  If
%       left blank, all genes will be shown.
% 
% OUTPUTS:
% 'DynGenesTable'
%       Table of dynamically expressed genes.  Table includes gene name,
%       fdr-corrected p-value, and order based on peak expression.
%

%% CODE:

% index genes and cells, order cells along trajectory position
cell_ind = find(~isnan(traj_pos));
traj_pos = traj_pos(cell_ind);
[~,cell_sort_ind] = sort(traj_pos);
X = full(X(:,cell_ind(cell_sort_ind)));

% filter out any genes below cells/counts threshold
min_max_flag = sum(X >= counts_thresh,2) >= min_cells;
X = X(min_max_flag,:);
gene_ind = find(min_max_flag);
nGenes = size(X,1);
 
% identify dynamically variable genes:
% (1) calculate p-value between two sliding windows w/ max and min average counts
nCells = length(cell_ind);
for j = 1:nGenes    
    % get sliding window average counts
    for k = 1:nCells-sliding_window+1
        wind(k,:) = k:k+sliding_window-1;
        tmp_X_avg(k) = mean(X(j,wind(k,:)));    
    end
    % get max and min windows
    [~,max_win_k]=max(tmp_X_avg);
    [~,min_win_k]=min(tmp_X_avg);
    % determine p-value in max vs min windows
    [~,p(j)] = ttest2(X(j,wind(max_win_k,:)),X(j,wind(min_win_k,:)));
end
% (2) calculate p-value again but with cell order randomized
rng(802); % sets the random number seed for reproducibility
rand_ind = randperm(nCells);
tmp_X_rand = X(:,rand_ind);
for j = 1:nGenes
    % get sliding window average counts
    for k = 1:nCells-sliding_window+1
        wind(k,:) = k:k+sliding_window-1;
        tmp_X_avg(k) = mean(tmp_X_rand(j,wind(k,:)));    
    end
    % get max and min windows
    [~,max_win_k]=max(tmp_X_avg);
    [~,min_win_k]=min(tmp_X_avg);
    % determine p-value in max vs min windows
    [~,p_rand(j)] = ttest2(tmp_X_rand(j,wind(max_win_k,:)),tmp_X_rand(j,wind(min_win_k,:)));
end
% (3) calculate fdr
for j = 1:nGenes
    fdr(j) = (sum(p_rand <= p(j))/nGenes);
    fdr_flag(j) = fdr(j) <= fdr_alpha;
end

% subset to only the dynamically variable genes
gene_ind = gene_ind(fdr_flag);
X = X(fdr_flag,:);

% zscore and smooth data
X_smooth = smoothdata(zscore(X')', 2, 'gaussian', sliding_window)';

% reorder genes by peak expression
[~, peak_cell_ind] = max(X_smooth);
[peak_cell_sorted, gene_ord] = sort(peak_cell_ind, 'ascend');
X_smooth = flipud(X_smooth(:,gene_ord));

% list the dynamically expressed genes, in order
dyn_genes = gene_names(gene_ind(gene_ord));

% subset which genes to show on the plot
if exist('gene_names_show', 'var')
    % first convert gene names to row indices
    for j = 1:length(gene_names_show)
        [~, gene_ind_show(j)] = ismember(gene_names_show{j}, gene_names);
    end
    cell_string_for_plot = repmat({''},1,length(gene_names));
    cell_string_for_plot(gene_ind_show) = gene_names(gene_ind_show);    
    dyn_genes_show = cell_string_for_plot(gene_ind(gene_ord));
else % show all genes
    dyn_genes_show = dyn_genes;
end

% generate table of dyn genes identified, p-values, and fdr flags
DynGenesTable = table(gene_names(gene_ind(gene_ord))', fdr(fdr_flag)', (1:length(dyn_genes))', peak_cell_sorted', 'VariableNames', {'Gene','FDR','Order', 'PeakCell'});

% make figure
figure('Position', [400 400 600 600])

% plot heatmap
ax1 = axes('position', [0.05, 0.04, 0.8, 0.9]);
imagesc(ax1, X_smooth(:,:)')
yTick_row_names = fliplr(dyn_genes_show);
set(gca,'YAxisLocation','right');
set(gca,'YTick',([1:length(yTick_row_names)]),'YTickLabel',yTick_row_names,'TickLength',[0 0])
set(gca,'XTickLabel','','TickLength',[0 0])
%colormap(ax1,(othercolor('BrBG10',1000)))
colormap(ax1,flipud(othercolor('RdGy10',1000)))
ax1.CLim = [-1 1];
%colorbar

% add pseudotime color bar
ax2 = axes('position', [0.05, 0.025, 0.8, 0.01]);
imagesc(ax2, sort(traj_pos)')
%colormap(ax2, othercolor('RdYlBu10',1000)); 
colormap(ax2, jet(1000));
%caxis(ax2,[0 0.9])
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0])
axis off

% link axes for zooming
linkaxes([ax1, ax2],'x');

