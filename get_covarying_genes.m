function gene_ind = get_covarying_genes(X_use, gene_ind, minimum_correlation)
% Usage: gene_ind = get_covarying_genes(X_use, gene_ind, minimum_correlation)
%
% Filters a gene list to include genes that are minimally correlated to at
% least one other gene.
%
% INPUTS:
% X_use                 Data matrix (rows=genes, columns=cells) for 
%                       constructing gene-gene correlations
% gene_ind              Row indices of genes to test
% minimum_correlation   Genes whose nearest neighbors are correlated below
%                       this threshold will be discarded
%
% OUTPUTS:
% gene_ind              Row indices of the filtered gene list
%

%% CODE:
warning('off')

nGenes = size(gene_ind,1);
X_ind = X_use(gene_ind,:);

% filter the top counts value for each gene.  avoids basing a correlation
% on a single outlier cell
[~,I] = max(X_ind,[],2);
for j = 1:nGenes
    X_ind(j,I(j)) = 0;
end

% compute gene-gene correlation matrix
gene_correlation_matrix = 1-pdist(X_ind, 'correlation');
gene_correlation_matrix = squareform(gene_correlation_matrix);

% only keep genes whose nearest neighbor is above threshold 
ind_keep = max(gene_correlation_matrix) > minimum_correlation;
gene_ind = gene_ind(ind_keep);

warning('on')

