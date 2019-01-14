function [table_CorrGenes, corr_gene_ind] = get_gene_neighbors(X, gene_names_all, gene_query, count_thresh, cell_thresh, ndisp)
% Usage: [table_CorrGenes, corr_gene_ind] = get_gene_neighbors(X, gene_names_all, gene_query, count_thresh, cell_thresh, ndisp)
%
% Given a single query gene, returns table of most highly correlated genes 
% based on Pearson distance, with associated P-values.
% 
% INPUTS:
% X        
%       Normalized matrix of single-cell expression counts (rows=genes, 
%       columns=cells).
%
% gene_names_all       
%       Names of all genes (cell array of strings), corresponding to 
%       rows of X.
%
% gene_query
%       User specified query gene (string).
% 
% counts_thresh & cell_thresh
%       Limit analysis to genes with at least counts_thresh number of 
%       counts in at least cell_thresh number of cells
%
% ndisp
%       Number of nearest gene neighbors to return.       
% 
% OUPUTS:
% table_CorrGenes
%       Table of height ndisp, summarizing the results of the neighbor
%       search.  Neighbor genes are rows.  Columns are: RHO (correlation),
%       PVAL, and adjusted PVAL (Benjamini & Hochberg 1995, fdr_bh).
%
% corr_gene_ind
%       Row indices of correlated genes identified.  
%

%% CODE:

% add path for p-value adjustment function fdr_bh
addpath('scTools/fdr_bh')

% extract the counts values for the query gene across all cells
query_gene_ind = strcmp(gene_names_all, gene_query);
query_gene_vals = X(query_gene_ind,:)';

% exclude genes without at least 'count_thresh' counts in 'cell_thresh cells'
% but still retain the query gene
gene_filt = sum(X>count_thresh,2)>cell_thresh | query_gene_ind';
X = X(gene_filt,:);
gene_names_all = gene_names_all(gene_filt);

% transpose of the original data such that rows=cells and columns=genes
Xtranspose = X';

% compute the correlation between selected gene & all other genes
[RHO,PVAL] = corr(Xtranspose,query_gene_vals,'rows','complete');
[~,~,~,PVAL_adj] = fdr_bh(PVAL);

% check how many genes to display
if ~exist('ndisp', 'var')
    ndisp = length(RHO);
end

% sort genes by correlation
RHO(isnan(RHO)) = -Inf;
[RHO_sorted,corr_gene_ind] = sort(RHO,'descend');
RHO_sorted = RHO_sorted(1:ndisp);
corr_gene_ind = corr_gene_ind(1:ndisp);
gene_names_all = gene_names_all(corr_gene_ind);

% generate results table
table_CorrGenes = table(RHO_sorted, PVAL(corr_gene_ind), PVAL_adj(corr_gene_ind),'VariableNames',{'RHO' 'PVAL' 'AdjPVAL'},'RowNames',(gene_names_all));


