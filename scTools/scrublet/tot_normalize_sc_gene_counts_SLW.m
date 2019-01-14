function [Xnm, counts, counts_orig, keep_gene_idx] = tot_normalize_sc_gene_counts_SLW(X, exclude_dominant_genes, max_frac)
% [Xnm, counts, counts_orig, keep_gene_idx] = tot_normalize_sc_gene_counts(X,true,0.1)

% INPUT
%
% X: raw counts


% OPTIONAL INPUT
%
% exclude_dominant_genes (boolean, default=false) and max_frac (float):
% If remove_dominant_genes is true, then any gene making up more than 
% max_frac of any cell's total counts is excluding from the total count 
% normalization for all cells.


% OUTPUT
%
% Xnm = matrix size(X) with normalized gene expression counts, such that
%           sum(Xnm(:,j)) = mean(counts) for all j.
%
% counts = vector size(X,2) with total counts per cell (corrected if
% excluding dominant genes)
%
% counts_orig = vector size(X,2) with total original counts per cell
% (including any excluded genes)
%
% keep_gene_idx = vector of size(X,2) with indices of the genes included in
% the normalization

counts = sum(X,1);% Total counts per cell
counts_orig = counts;
keep_gene_idx=true(size(X,1),1);

if(exist('exclude_dominant_genes','var'))
    if exclude_dominant_genes
        for iCell=1:size(X,2)
            gene_fracs=X(:,iCell)/counts(iCell);
            keep_gene_idx(gene_fracs>max_frac) = false;
        end
        counts=sum(X(keep_gene_idx,:),1);
    end
end

w = mean(counts)./counts; % normalization factor
tmp = repmat(w,size(X,1),1);
Xnm = X.*tmp;


