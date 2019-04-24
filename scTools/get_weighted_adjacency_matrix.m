function [Aw_sym, Aw_non_sym] = get_weighted_adjacency_matrix(G)
% Usage: [Aw_sym, Aw_non_sym] = get_weighted_adjacency_matrix(G)
%
% Converts a Matlab graph object into a symmetric adjacency matrix with
% edge weights.  
%
%% CODE: 

nNodes = numnodes(G);
[s,t] = findedge(G);
Aw_non_sym = sparse(s, t, G.Edges.Weights, nNodes, nNodes);
Aw_sym = Aw_non_sym + Aw_non_sym';

