function [X_norm, tot_counts] = get_normalized_counts(X)
% Usage: [X_norm, tot_counts] = get_normalized_counts(X)
%
% Perform total counts normalization.
%
% INPUTS
% X            Raw UMI-filtered counts matrix (rows=genes, columns=cells)
%
% OUTPUTS
% X_norm       Normalized counts matrix
% tot_counts   Array of original UMI count totals for each cell (column) in
%              X

%% CODE:

tot_counts = sum(X,1);
w = mean(tot_counts)./tot_counts; 
tmp = repmat(w,size(X,1),1);
X_norm = X.*tmp;


