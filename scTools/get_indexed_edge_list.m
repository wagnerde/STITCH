function edge_list = get_indexed_edge_list(G)
% Usage: edge_list = get_indexed_edge_list(G)
%
% Generates an edge list from graph object G.  Endnodes are renumbered 
% from 0.  Duplicate (mutual) edges are removed.
%
%% CODE:

% get indexed source and target nodes for each edge
full_adj = full(adjacency(G));
[source_nodes, target_nodes] = find(full_adj);
edge_list = [source_nodes, target_nodes] - 1; % index to zero  

% mutual edges are now duplicate rows, non-mutual edges will have only a single row
[~,ia,ic] = unique(sort(edge_list(:,[1,2]), 2 ), 'rows'); 

% keep only the first edge from each duplicate pair
nOccurences = histc(ic,unique(ic));
is_duplicate_edge = (nOccurences == 2); 
keep_ind = ia(is_duplicate_edge);  
edge_list = fliplr(edge_list(keep_ind,:)); 