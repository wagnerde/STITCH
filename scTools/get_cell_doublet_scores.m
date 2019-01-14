function doublet_scores = get_cell_doublet_scores(X_norm, tot_counts, gene_ind, doub_frac, k_neighbors)
% Usage: doublet_scores = get_cell_doublet_scores(X, tot_counts, gene_ind, doub_frac, k_neighbors)
%
% Executes a Matlab version of "Scrublet", a cell doublet filtering
% algorithm by Sam Wolock, see also:
% https://www.biorxiv.org/content/early/2018/07/09/357368
% https://github.com/AllonKleinLab/scrublet
%
% Scrublet includes its own private functions for total counts normalization 
% and variable gene selection. This implementation relies on pre-calculating
% these objects with scTools, for consistency.
%

%% CODE:

% add functions to path:
addpath('scTools/scrublet')

% get doublet scores
doublet_scores = doublet_detector(X_norm, 'counts', tot_counts, 'genes_use', gene_ind, 'doub_frac', doub_frac, 'k_neighbors', k_neighbors);



