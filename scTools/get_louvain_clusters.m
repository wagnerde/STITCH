function [ClustID, ClustID_filt, Quality] = get_louvain_clusters(G, gamma)
% Usage: [ClustID, ClustID_filt, Quality] = get_louvain_clusters(G, gamma)
%
% Performs Louvain-based community detection using 'genlouvain'
%
% INPUTS:
% G         Graph object
% gamma     Resolution (default = 1). Values >1 gives more clusters.
%           Values <1 give fewer clusters.
% 
% OUTPUTS:
% ClustID         Community/Cluster assignments for each cell
% ClustID_filt    Small clusters assigned to NaN
% Quality         Quality score for the partition of the network
%
%% CODE:

% add function path
addpath('scTools/GenLouvain-master')

% compute jaccard weights for each edge
A = get_jaccard_edge_weights(G);

% compute the modularity/quality matrix B, perform community detection
k = full(sum(A));
twom = sum(k);
B = @(v) A(:,v) - gamma*k'*k(v)/twom;
[ClustID, Quality] = genlouvain(B,[],0);
Quality = Quality/twom;

% ignore clusters with < 10 cells
ClustID_filt = ClustID;
[clust_counts, ~] = count_unique(ClustID_filt);
low_count_flag = find(clust_counts < 10);
ClustID_filt(ismember(ClustID_filt, low_count_flag)) = NaN;

% clean up cluster assignments
ClustID_filt = grp2idx(ClustID_filt);
ClustID = grp2idx(ClustID);

