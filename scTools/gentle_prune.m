function G = gentle_prune(G, edges)
% Usage: G = gentle_prune(G, edges)
% 
% Accepts a graph and a list of candidate edges to prune
% Sort edges by increasing weight. Attempt to prune all edges, escaping if
% doing so would cut the graph
%
%% SETTINGS:
verbose = false;

%% CODE:

% sort edges 
[~ , ind] = sort(G.Edges.Weights(edges,1));
edges_sorted = edges(ind);
edge_sorted_s = G.Edges.EndNodes(edges_sorted,1);
edge_sorted_t = G.Edges.EndNodes(edges_sorted,2);

% number of starting connected components
starting_cc = length(unique(conncomp(G)));

% cut edges if doing so does not increase the # of connected components
for k = 1:length(edges_sorted)
    G_cut_test = rmedge(G, edge_sorted_s(k), edge_sorted_t(k));    
    if length(unique(conncomp(G_cut_test))) == starting_cc
        G = rmedge(G, edge_sorted_s(k), edge_sorted_t(k));  
    elseif verbose
        disp(['Avoided severing branch at edge ' num2str(k)]);        
    end
end