function node_IDs = stitch_celldata_to_nodeID(DataSet, G, celldata_var)
% Usage: node_IDs = stitch_celldata_to_nodeID(DataSet, G, celldata_var)
%
% Indexes the 'celldata' annotations table of DataSet to STITCH graph G
% 
% INPUTS:
% DataSet           STITCH dataSet object
% G                 STITCH graph object
% cell_data_var     String. Must match one of the celldata column/variable
%                   names. Specifies the set of celldata annotations which
%                   will be used for indexing.
% OUTPUTS:
% node_IDs          String array of annotations for each node in G.

%% CODE: 

node_IDs = [];
celldata = {DataSet.celldata};
for j = length(celldata):-1:1
    % determine which column of celldata to index
    celldata_var_ind = strcmp(celldata{1}.Properties.VariableNames, celldata_var);
    % concatenate node_IDs for each DataSet entry in reverse time order
    node_IDs = [node_IDs; string(table2cell(celldata{j}(:,celldata_var_ind)))]; 
end

% index node_IDs to graph G
node_IDs = index_to_graph(G, node_IDs);
