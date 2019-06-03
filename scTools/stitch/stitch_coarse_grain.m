function [G_cg, G_cg_scaff] = stitch_coarse_grain(G, node_IDs, varargin)
% Usage: [G_cg, G_cg_scaff] = stitch_coarse_grain(G, node_IDs, varargin)
%  
% Generates a coarse-grained ("cg") kNN graph by collapsing nodes with the
% same group IDs.
%
% INPUTS:
% G                STITCH graph object (nodes=single cells)
% node_IDs         String array of identifiers (e.g. cluster names) for 
%                  each node in G
% 
% Optional input name/value pairs: 
% 'edge_weight_norm' 
%                  Acceptable inputs are 1 or 2 (default=1). Sets the 
%                  strategy for normalizing edge weights in the 
%                  coarse-grained graph. In both cases, the number of 
%                  original shared single cell edges is scaled to determine
%                  weights of the resultingcoarse-grained edge.
%                  1: # of shared edges is scaled to the total number of
%                  outgoing edges across the two participating node_IDs
%                  (e.g. the Jaccard Index)
%                  2: # of shared edges is scaled to the number of outgoing 
%                  edges from the less abundant of the two node_IDs
%                   
% 'edge_weight_thresh'
%                  Minimum normalized edge weight (default=0.01).  
%                  Coarse-grained edges with normalized weights below this 
%                  value are filtered.
%
% 'edge_weight_histogram' 
%                  Boolean (default=false).  Plot a histogram of normalized
%                  edge weights with 'edge_weight_thresh' filter indicated.
% 
% OUTPUTS:
% G_cg             Graph object corresponding to the coarse-grained graph.
% G_cg_scaff       A subgraph of G_cg containing only scaffold edges (e.g.
%                  a spanning tree w/ most heavily weighted edges
%                  connecting each coarse-grained node to the previous 
%                  timepoint.
%

%% PARAMETER SETTINGS
% Set defaults
def.edge_weight_thresh = 0.01;
def.edge_weight_norm = 1;
def.edge_weight_histogram = false;

% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'stitch_coarse_grain';
parserObj.StructExpand = false; 
parserObj.addOptional('edge_weight_thresh',def.edge_weight_thresh);
parserObj.addOptional('edge_weight_norm',def.edge_weight_norm);
parserObj.addOptional('edge_weight_histogram',def.edge_weight_histogram);

% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% CODE:

% if necessary, index node_IDs to G
if length(G.Nodes.Name) ~= length(node_IDs)
    node_IDs = index_to_graph(G, node_IDs);
end

% flag any unannotated nodes
unannotated_nodes = ismissing(node_IDs);
node_IDs(unannotated_nodes) = "Unassigned";

% convert node_IDs to numeric
node_IDs_orig = node_IDs;
node_IDs = grp2idx(node_IDs);

% get all unique node_IDs and their names contained in the graph
unique_nodeIDs_in_graph = unique(node_IDs, 'stable');
unique_node_names_in_graph = unique(node_IDs_orig, 'stable');
unannotated_node_ID = find(strcmp(unique_node_names_in_graph,'Unassigned'));
unique_node_names_in_graph(unannotated_node_ID) = [];
unique_nodeIDs_in_graph(unannotated_node_ID) = [];
nCGNodes = max(unique_nodeIDs_in_graph);

% get timepoints for each node_ID
for j=unique_nodeIDs_in_graph'
    nodeID_times(j) = mean(G.Nodes.OriginalDataSet(node_IDs==j));
end
nodeID_times_adj = nodeID_times;
nodeID_times_adj(unannotated_node_ID)=[];

% build a coarse-grained edge table
% nodes = unique node_IDs 
% edges = weighted by # of original edges connecting those two node_IDs
CG_EdgeTable = [];
j1 = 1;
for j = unique_nodeIDs_in_graph' % ...for each new "CG" node...
    % get CG node size (number of original cells with this node_ID)
    CG_node_sizes(j1,1) = sum(node_IDs == j);   
    % get original single-cell node indices of this particular node_ID
    node_indices = find(node_IDs == j);   
    % get all neighbors (single-cell indices) to this node_ID
    neighborhood = [];
    for k = 1:length(node_indices)
        neighborhood = [neighborhood; neighbors(G, node_indices(k))];
    end    
    % convert old single-cell node indices to node_IDs
    nodeID_neighborhood = node_IDs(neighborhood);   
    % remove internal/self edges
    nodeID_neighborhood = nodeID_neighborhood(nodeID_neighborhood ~= j);    
    % remove any edges to the "unannotated" / NaN node_ID
    if ~isempty(unannotated_node_ID)
        nodeID_neighborhood = nodeID_neighborhood(nodeID_neighborhood ~= unannotated_node_ID);   
    end
    % count all unique outgoing (non-self) edges
    non_self_targets = unique(nodeID_neighborhood);
    nOriginalEdges = histc(nodeID_neighborhood, non_self_targets);
    nOriginalEdges_for_this_NodeID(j1) = sum(nOriginalEdges);    
    % flag this edge depending on whether it connects within or between timepoints
    nCGEdges = length(non_self_targets);
    self = repmat(j,1,nCGEdges)';
    is_within_timepoint = nodeID_times(self) == nodeID_times(non_self_targets);
    is_between_timepoint = nodeID_times(self) ~= nodeID_times(non_self_targets);    
    % append EdgeTable for this clustID
    collapsed_edge_table_next = table([self, non_self_targets], zeros(length(nOriginalEdges),1), nOriginalEdges, zeros(length(nOriginalEdges),1), is_within_timepoint', is_between_timepoint', 'VariableNames', {'EndNodes', 'Weights', 'nOriginalEdges', 'Denom', 'WithinTimePoint', 'BetweenTimePoint'});    
    CG_EdgeTable = [CG_EdgeTable; collapsed_edge_table_next];
    j1=j1+1;
end

% build a coarse-grained node table
CG_NodeTable = [];
CG_NodeTable = table(unique_nodeIDs_in_graph, unique_node_names_in_graph, CG_node_sizes, nodeID_times_adj', nodeID_times_adj', nOriginalEdges_for_this_NodeID', 'VariableNames', {'Name', 'NodeLabels', 'nCells', 'Timepoint', 'OriginalDataSet', 'nOriginalOutgoingEdges'});

% calculate fraction of cells from each timepoint w/ each node_ID
unique_nodeID_times = unique(nodeID_times_adj);
for j = 1:length(unique_nodeID_times)    
    nodeID_in_this_timepoint = nodeID_times(unique_nodeIDs_in_graph) == unique_nodeID_times(j);
    this_timepoint_nCells(j) = max(sum(CG_NodeTable.nCells(nodeID_in_this_timepoint)),1);
end
for j = 1:height(CG_NodeTable)
    CG_NodeTable.fCells(j,1) = CG_NodeTable.nCells(j,1) / this_timepoint_nCells(CG_NodeTable.Timepoint(j,1));
end

% remove duplicate edges
CG_EdgeTable = filter_duplicate_edges(CG_EdgeTable);

% give each edge a unique identifier
CG_EdgeTable.UniqueID = (1:height(CG_EdgeTable))';

% normalize edge weights according to 'settings.edge_weight_norm'
CGNodeNames = CG_NodeTable.Name;
[~,Node1_row_ind] = ismember(CG_EdgeTable.EndNodes(:,1),CGNodeNames);
[~,Node2_row_ind] = ismember(CG_EdgeTable.EndNodes(:,2),CGNodeNames);
if settings.edge_weight_norm == 2 
    % denominator = # of outgoing edges in the smaller of the two clusters
    CG_EdgeTable.Denom = min(([CG_NodeTable.nOriginalOutgoingEdges(Node1_row_ind) CG_NodeTable.nOriginalOutgoingEdges(Node2_row_ind)]),[],2);
elseif settings.edge_weight_norm == 1 
    % denominator = sum total # of outgoing edges in both clusters
    CG_EdgeTable.Denom = sum(([CG_NodeTable.nOriginalOutgoingEdges(Node1_row_ind) CG_NodeTable.nOriginalOutgoingEdges(Node2_row_ind)]),2);
end
CG_EdgeTable.Weights = CG_EdgeTable.nOriginalEdges ./ CG_EdgeTable.Denom;

% flag edges according to latest timepoint of either participating node
CG_EdgeTable.EdgeTimes = max(nodeID_times(CG_EdgeTable.EndNodes(:,1)), nodeID_times(CG_EdgeTable.EndNodes(:,2)))';

% generate coarse-grained graph object
CG_NodeTable.Name = cellstr(num2str(CG_NodeTable.Name));
CG_EdgeTable.EndNodes = [cellstr(num2str(CG_EdgeTable.EndNodes(:,1))) cellstr(num2str(CG_EdgeTable.EndNodes(:,2)))];
G_cg = graph(CG_EdgeTable, CG_NodeTable);

% convert NodeIDs to numerical array
G_cg.Edges.NodeIDs = [str2num(cell2mat(G_cg.Edges.EndNodes(:,1))),str2num(cell2mat(G_cg.Edges.EndNodes(:,2)))];
G_cg.Edges.NodeIDs = sort(G_cg.Edges.NodeIDs(:,[1,2]),2);

% plot edge weight histogram
if settings.edge_weight_histogram
    figure; 
    hist((G_cg.Edges.Weights),1000);   
    title('Edge Weight Histogram')
    xlabel('Edge Weights (log10)')
    ylabel('# of Edges')
    y_limits = ylim; 
    line([(settings.edge_weight_thresh) (settings.edge_weight_thresh)],[y_limits(1) y_limits(2)],'color','r','linewidth',1)
end

% sort edges starting from last timepoint and work backwards to the first timepoint
unique_edge_times = sort(unique(G_cg.Edges.EdgeTimes), 'descend');

% flag the weakest cg graph edges according to 'settings.edge_weight_thresh'
G_weak_rmvd = G_cg;
for j=1:length(unique_edge_times)
    weak_edges_in_this_time = (G_weak_rmvd.Edges.EdgeTimes == unique_edge_times(j)) & (G_weak_rmvd.Edges.Weights < settings.edge_weight_thresh);    
    G_weak_rmvd = rmedge(G_weak_rmvd, find(weak_edges_in_this_time));
end
G_cg.Edges.WeakRmvd = ismember(G_cg.Edges.UniqueID, G_weak_rmvd.Edges.UniqueID);

% generate a "scaffold" subgraph (a spanning tree w/ strongest edges 
% connecting to previous timepoint)
G_cg_scaff = G_cg;
for j=1:length(unique_edge_times)
    weak_edges_in_this_time = (G_cg_scaff.Edges.EdgeTimes == unique_edge_times(j)) & (G_cg_scaff.Edges.Weights < Inf);    
    G_cg_scaff = gentle_prune(G_cg_scaff, find(weak_edges_in_this_time));
end
G_cg.Edges.Scaffold = ismember(G_cg.Edges.UniqueID, G_cg_scaff.Edges.UniqueID);

% remove all edges below threshold
G_cg = rmedge(G_cg, find(~G_cg.Edges.Scaffold & ~G_cg.Edges.WeakRmvd));

% isolate the 'giant' component (connected component #1)
G_cg = filter_graph_giant_component(G_cg);
G_cg_scaff = filter_graph_giant_component(G_cg_scaff);


