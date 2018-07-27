function stitch_filtered = stitch_filter(stitch_merged, settings)
% Usage: stitch_filtered = stitch_filter(stitch_merged, settings)
% 
% Processes a set of merged node & edge tables. Applies a series of local 
% and gloval edge distance filters, enforces mutual edges, and returns the 
% largest graph connected component.
% 
% INPUTS:
% stitch_merged         A structure array containing a single merged node
%                       table and a single merged edge table
% settings              Structure containing STITCH settings
% 
% OUTPUT:
% stitch_filtered       A structure array containing the final filtered
%                       Matlab graph object and adjacency matrix
% 

%% CODE:

% load nodes and edges
EdgeTable = stitch_merged.EdgeTable;
NodeTable = stitch_merged.NodeTable;
disp('Filtering graph...');
disp(['# starting nodes: ', num2str(height(NodeTable))]);
disp(['# starting edges: ', num2str(round(height(EdgeTable)/1000)), 'K']);

% (1) filter based on hard thresholds for global and local neighbor distances
graph_edge_hard_thresh_fail = (EdgeTable.D_local > settings.graph_max_D_local | EdgeTable.D_global > settings.graph_max_D_global) & EdgeTable.InternalEdge;    
link_edge_hard_thresh_fail = (EdgeTable.D_local > settings.link_max_D_local | EdgeTable.D_global > settings.link_max_D_global) & EdgeTable.LinkEdge;
edge_hard_thresh_fail = graph_edge_hard_thresh_fail | link_edge_hard_thresh_fail;
EdgeTable(edge_hard_thresh_fail,:) = [];  
disp(['# nodes after distance filters: ', num2str(height(NodeTable))]);
disp(['# edges after distance filters: ', num2str(round(height(EdgeTable)/1000)), 'K']);

% (2) retain only the basis-outgoing edges that had a mutual incoming edge
if settings.require_mutual_edges
    EdgeTable = filter_nonmutual_edges(EdgeTable);
    EdgeTable(EdgeTable.BasisOutgoing == 0,:) = [];
    disp(['# nodes after mutual edge filter: ', num2str(height(NodeTable))]);
    disp(['# edges after mutual edge filter: ', num2str(round(height(EdgeTable)/1000)), 'K']);
else
    disp('Considering all outgoing edges')
    EdgeTable(EdgeTable.BasisOutgoing == 0,:) = [];
    disp(['# nodes remaining: ', num2str(height(NodeTable))]);
    disp(['# edges remaining: ', num2str(round(height(EdgeTable)/1000)), 'K']);
end

% (3) restrict graph to k_final outgoing edges per node
% For each node, consider all graph and link edges built in that node's 
% PCA basis space. Rank by distance, keep only the 'shortest' k_final edges. 
[~, edges_by_distance] = sort(EdgeTable.D_orig);
all_source_NodeID = EdgeTable.EndNodes(:,1);
source_NodeIDs_by_distance = all_source_NodeID(edges_by_distance);
all_uniqueID = EdgeTable.UniqueEdgeID;
uniqueID_by_distance = all_uniqueID(edges_by_distance);
% loop depending on how many edges we want to retain per node
UniqueEdgeID_to_keep = [];
for k = 1:settings.k_final
    % use the unique function to pull the indices for the first occurrence 
    % of each NodeID - these are the shortest edges
    [~, ia, ~] = unique(source_NodeIDs_by_distance);
    % append the uniqueID of edges we want to keep
    UniqueEdgeID_to_keep = [UniqueEdgeID_to_keep; uniqueID_by_distance(ia)];       
    % remove these edges from the by_distance lists for the next loop iteration
    source_NodeIDs_by_distance(ia) = [];   
    uniqueID_by_distance(ia) = [];   
end
UniqueEdgeID_to_keep_original_ind = ismember(EdgeTable.UniqueEdgeID, UniqueEdgeID_to_keep);
EdgeTable = EdgeTable(UniqueEdgeID_to_keep_original_ind,:);
EdgeTable = filter_duplicate_edges(EdgeTable);    
disp(['# nodes after k_final filter: ', num2str(height(NodeTable))]);
disp(['# edges after k_final filter: ', num2str(round(height(EdgeTable)/1000)), 'K']);

% (4) isolate the 'giant' component (connected component #1)
G = graph(EdgeTable, NodeTable); 
G = filter_graph_giant_component(G);
disp(['# nodes in giant: ', num2str(numnodes(G))]);
disp(['# edges in giant: ', num2str(round(numedges(G)/1000)), 'K']);

% calculate adjacency matrix
A = adjacency(G);

% export
stitch_filtered.G = G;
stitch_filtered.A = A;

disp('Done filtering');

