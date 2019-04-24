function gene_ind = get_correlated_genes(gene_ind, X, corr_thresh, nIter)
% Usage: gene_ind = get_correlated_genes(gene_ind, X, corr_thresh, nIter)
%
% Expands a list of genes by identifying additional genes with highly 
% correlated single-cell expression patterns.
% 
% INPUTS:
% gene_ind      Row indices of genes to seed the search
% X             Data matrix (rows=genes, columns=cells) for constructing
%               single-cell gene-gene correlations
% corr_thresh   Minimum correlation for adding new genes
% nIter         Number of rounds of list expansion
%
% OUTPUTS:
% gene_ind      Row indices of the expanded gene list
% 

%% CODE

for j = 1:nIter
        
    c = corr(X(gene_ind,:)',X');

    for k = 1:length(gene_ind)
        gene_ind = union(gene_ind,find((c(k,:) >= corr_thresh))); 
    end

end
