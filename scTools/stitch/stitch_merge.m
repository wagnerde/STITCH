function stitch_merged = stitch_merge(stitch_knn_graphs, stitch_links)
% Usage: stitch_merged = stitch_merge(stitch_knn_graphs, stitch_links)
% 
% Merges together a series of node and edge tables. All NodeIDs are 
% updated with unique names in reverse chronological order  
%
% INPUTS:
% stitch_knn_graphs   A structure array of all individual within-
%                     timepoint node & edge tables
% stitch_links        A structure array of all link node & edge tables
%
% OUTPUT:
% stitch_merged       A structure array containing a single merged node
%                     table and a single merged edge table
%

%% CODE: 

% determine number of graph and link tables
nGraphs = length(stitch_knn_graphs);
nLinks = length(stitch_links);

% generate tally and base_counter arrays
% tally is the number of cells in each timepoint
% base_counter is the cumulative sum over timepoints
build_order = cell2mat({stitch_knn_graphs.DataSet_ind});
tally = cell2mat({stitch_knn_graphs.nNodes});
tmp = 0;
for j = 1:length(tally)
    base_counter(j) = tmp;
    tmp = tmp + tally(j);
end

% stitch node and edge tables from each kNN graph
% increment node indices by base_counter 
EdgeTable_stitched = stitch_knn_graphs(1).EdgeTable;
NodeTable_stitched = stitch_knn_graphs(1).NodeTable;

for j = 2:nGraphs

    EdgeTable_next = stitch_knn_graphs(j).EdgeTable;
    EdgeTable_next.EndNodes = EdgeTable_next.EndNodes + base_counter(j);
    EdgeTable_stitched = [EdgeTable_stitched; EdgeTable_next];

    NodeTable_next = stitch_knn_graphs(j).NodeTable;
    NodeTable_next.Name = NodeTable_next.Name + base_counter(j);
    NodeTable_stitched = [NodeTable_stitched; NodeTable_next];

end

% stitch link tables
% each set of link tables will be from arbitrary basis and link datasets
% node indices need to be incremented accordingly
for j = 1:nLinks
    
    basis_ind = build_order(stitch_links(j).basis_ind);
    link_ind = build_order(stitch_links(j).link_ind);
    
    basis_tally = tally(basis_ind);
    link_tally = tally(link_ind);
    basis_base = base_counter(basis_ind);
    link_base = base_counter(link_ind);
    
    EdgeTable_next = stitch_links(j).EdgeTable;
    
    % check to make sure the node total in the Edge table does not exceed the size of the basis and link tables
    if max(EdgeTable_next.EndNodes(:)) > basis_tally + link_tally
        disp(['Error! Node sum mismatch on link ' num2str(j)])
    end
    
    % update node indices
    % basis nodes will ALWAYS be listed first in the link EdgeTable (i.e. they will initially have lower NodeIDs)
    basis_nodes = EdgeTable_next.EndNodes <= basis_tally;
    link_nodes = EdgeTable_next.EndNodes > basis_tally;    
    % basis NodeIDs are currently 1:basis_tally -> just add base_counter
    EdgeTable_next.EndNodes(basis_nodes(:,1),1) = EdgeTable_next.EndNodes(basis_nodes(:,1),1) + basis_base;
    EdgeTable_next.EndNodes(basis_nodes(:,2),2) = EdgeTable_next.EndNodes(basis_nodes(:,2),2) + basis_base;
    % link NodeIDs are currently (1:link_tally)+basis_tally -> first subtract basis_tally then just add link_base
    EdgeTable_next.EndNodes(link_nodes(:,1),1) = EdgeTable_next.EndNodes(link_nodes(:,1),1) - basis_tally + link_base;
    EdgeTable_next.EndNodes(link_nodes(:,2),2) = EdgeTable_next.EndNodes(link_nodes(:,2),2) - basis_tally + link_base;
    
    % merge the current EdgeTable
    EdgeTable_stitched = [EdgeTable_stitched; EdgeTable_next];
    
end

% format node names as array of cell strings
NodeTable_stitched.Name = cellstr(num2str(NodeTable_stitched.Name));

% give each edge a unique ID
EdgeTable_stitched.UniqueEdgeID = (1:height(EdgeTable_stitched))';

% outputs
stitch_merged.NodeTable = (NodeTable_stitched);
stitch_merged.EdgeTable = (EdgeTable_stitched);
stitch_merged.tally = tally;
stitch_merged.base_counter = base_counter;

disp('Done merging');




