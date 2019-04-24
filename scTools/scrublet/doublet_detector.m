function [doub_score_obs, doub_score_full, doub_labels, PCdat] = doublet_detector(E, varargin)

%% OUTPUTS:
% doub_score_obs: vector of doublet scores (between 0 and 1) for each input
% cell. This is the main output; the rest are mostly used for testing and
% troubleshooting.
%    size: (number input cells, 1)

% doub_score_full: vector of doublet scores (between 0 and 1) for both
% input cells and simulated doublets
%    size: (number input cells + number simulated doublets, 1)

% doub_labels: vector of 0's and 1's labeling each entry of doub_score_full
% (0 = input cell; 1 = simulated doublet)
%    size: (number input cells + number simulated doublets, 1)

% PCdat: PCA coordinates for input cells and simulated doublets
%    size: (number input cells + number simulated doublets, number of PCs)

%% INPUTS

% E: Matrix of raw transcript counts (preferably unnormalized)
%    Rows are genes, columns are cells


%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional arguments:

% doub_frac: number of doublets to simulate, as a fraction of the total
% number of observed cells. (default = 3)

% k_neighbors: number of neighbors for kNN graph. This will automatically
% be scaled by doub_frac to prevent doublet saturation. A rule of thumb is
% to set k_neighbors=sqrt(number of input cells), though the optimal value
% likely depends on the true doublet rate and the cell population
% structure. (default = 50)

% distance_metric: distance metric for building kNN graph (default =
% 'euclidean')

% counts: total transcript counts for each input cell; most useful if
% you've already normalized the expression matrix E or are using the
% optional precomputed_pca input. (default = [])

% precomputed_pca: PCA coordinates for each input cell (cells are rows, PCs
% are columns). It's often best to run the doublet detector using the same
% PCA coordinates you've used for visualization purposes (e.g. as the input
% to t-SNE or SPRING). If you provide precomputed_pca, then none of the
% remaining optional arguments will be used, and the counts matrix E will
% only be used for calculating transcript count totals. (default = [])

% num_pc: number of principal components to use as input building the kNN
% graph (default = 50)

% exclude_dominant_gene_frac: when normalizing the expression matrix using
% total counts normalization, it is sometimes beneficial to exclude genes
% that make up a very large fraction of a cell's total counts. This
% argument sets this maximum fraction. (default = 1.0)

% genes_use: array of gene indices (rows of E) to use as input to PCA. This
% overrides all remaining options. (default = [])

% min_counts and min_cells: For filtering genes based on expression levels.
% To be included, genes must be expressed in at least min_counts copies in
% at least min_cells cells. (for both options, default = 3)

% vscore_percentile: Float between 0 and 1. For filtering genes based on
% variability. V-score is a measure of above-Poisson noise. To be included,
% genes must have a V-score in the top vscore_percentile percentile.
% (default = 0.85)


%% Setup
% Set defaults
def.doub_frac = 3;
def.k_neighbors = 50;
def.distance_metric = 'euclidean';
def.counts = [];
def.precomputed_pca = [];
def.num_pc = 50;
def.exclude_dominant_gene_frac = 1;
def.genes_use = [];
def.min_counts = 3;
def.min_cells = 3;
def.vscore_percentile = 0.85;


% Create parser object:
parserObj = inputParser;
parserObj.FunctionName = 'doublet_detector';
parserObj.StructExpand = false;
parserObj.addOptional('counts',def.counts);
parserObj.addOptional('precomputed_pca',def.precomputed_pca);
parserObj.addOptional('exclude_dominant_gene_frac',def.exclude_dominant_gene_frac);
parserObj.addOptional('genes_use',def.genes_use);
parserObj.addOptional('min_counts',def.min_counts);
parserObj.addOptional('min_cells',def.min_cells);
parserObj.addOptional('vscore_percentile',def.vscore_percentile);
parserObj.addOptional('num_pc',def.num_pc);
parserObj.addOptional('doub_frac',def.doub_frac);
parserObj.addOptional('k_neighbors',def.k_neighbors);
parserObj.addOptional('distance_metric',def.distance_metric);



% Parse input options:
parserObj.parse(varargin{:});
opt = parserObj.Results;

%%

if isempty(opt.counts)
    counts = sum(E, 1);
else
    counts = opt.counts;
end

if isempty(opt.precomputed_pca)

% DEW: supply a pre-normalized data matrix
%
%     disp('Total count normalizing')
%     if opt.exclude_dominant_gene_frac < 1
%         E = tot_normalize_sc_gene_counts_SLW(E, true, opt.exclude_dominant_gene_frac);
%     else
%         E = tot_normalize_sc_gene_counts_SLW(E, false);
%     end
%     
    if isempty(opt.genes_use)
        disp('Finding highly variable genes')
        v_scores = get_single_cell_Vscores(E, counts, 'show_plot',false);
        gene_filter = find((sum(E>=opt.min_counts, 2) >= opt.min_cells) & (v_scores>quantile(v_scores,opt.vscore_percentile)));
    else
        gene_filter = opt.genes_use;
    end
    
    disp(['Using ' num2str(length(gene_filter)) ' genes for PCA'])
    [~, PCdat] = pca(zscore(E(gene_filter,:), [], 2)', 'NumComponents', opt.num_pc);
    
    disp('Simulating doublets')
    [PCdat, doub_labels, ~] = simulate_doublets_from_pca(PCdat, counts, opt.doub_frac);
    
else
    PCdat = opt.precomputed_pca;
    disp('Simulating doublets')
    [PCdat, doub_labels, ~] = simulate_doublets_from_pca(PCdat, counts, opt.doub_frac);
end

n_obs = sum(doub_labels == 0);
n_sim = sum(doub_labels == 1);

k_detect = round(opt.k_neighbors * (1 + n_sim / n_obs));

disp(['Running KNN classifier with k = ' num2str(k_detect)])

% find k nearest neighbors for each cell (include both observed cells
% and simulated doublets in the search)
knn = knnsearch(PCdat,PCdat,'k',k_detect + 1,'distance',opt.distance_metric);

% remove first column (closest neighbor, which is the cell itself)
knn = knn(:,2:end);

n_doub_neigh = sum(doub_labels(knn) == 1, 2);
n_sing_neigh = sum(doub_labels(knn) == 0, 2);

doub_score_full = n_doub_neigh ./ (n_doub_neigh + n_sing_neigh * n_sim / n_obs);
doub_score_obs = doub_score_full(doub_labels == 0);

disp('Done.')


%%
function [PCdoub, doub_labels, pair_ix] = simulate_doublets_from_pca(PCdat, total_counts, doublet_frac)
    n_observed = size(PCdat, 1);
    n_doublet = round(n_observed * doublet_frac);
    
    pair_ix = randi(n_observed, n_doublet, 2);
    pair_tots = total_counts(pair_ix);
    pair_fracs = bsxfun(@rdivide, pair_tots, sum(pair_tots, 2));
    
    PCdoub = bsxfun(@times, PCdat(pair_ix(:, 1), :), pair_fracs(:, 1)) + bsxfun(@times, PCdat(pair_ix(:, 2), :), pair_fracs(:, 2));
    PCdoub = vertcat(PCdat, PCdoub);
    doub_labels = vertcat(zeros(n_observed, 1), ones(n_doublet, 1));

end

end