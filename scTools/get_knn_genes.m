function kNN_struct = get_knn_genes(X_norm, gene_ind, gene_names_all, k_use, distance_metric)
% Usage: kNN_struct = get_knn_genes(X_norm, gene_ind, gene_names_all, k_use, distance_metric)
%
% Identify k-nearest neighbor GENES from a normalized single-cell counts 
% matrix with knnsearch, then return node/edge tables. 
%
% INPUTS:
% 'X_norm'      
%         A total-counts normalized single-cell UMI counts matrix
%         rows=transcripts, columns=cells
%
% 'gene_ind'    
%         Row indices for genes of interest (e.g. highly variable genes)  
%
%
% 'k_use' 
%         Number of nearest neighbors to retain for each cell.  
%
% 'distance_metric' 
%         Distance metric to supply to knnsearch (e.g. 'euclidean').
% 
% OUTPUTS:
% kNN_struct 
%         Structure array containing node tables and edge tables
%
%

%% CODE

% subset the data matrix
X_ind = X_norm(gene_ind,:);
nNodes = size(X_ind,1);

% search for neighbors
k_use = k_use + 1;
[IDX_node2, D] = knnsearch(X_ind, X_ind, 'K', k_use, 'Distance', distance_metric);
IDX_node1 = repmat((1:(nNodes))',1,k_use);

% convert matrices to 1d arrays; omit column1 (self-edges)
IDX_node2 = mat2array(IDX_node2(:,2:end));
IDX_node1 = mat2array(IDX_node1(:,2:end));
D = mat2array(D(:,2:end));

% convert nodeIDs back to their original gene names
NodeIDs = gene_names_all(gene_ind((1:nNodes)'));

% assemble edge & node tables
EdgeTable = table([IDX_node1, IDX_node2], D, 'VariableNames', {'EndNodes', 'Weights'});
NodeTable = table(NodeIDs, 'VariableNames', {'Name'});
%NodeTable.Name = cellstr((NodeTable.Name));

% outputs
kNN_struct.nEdges = length(IDX_node2);
kNN_struct.nNodes = nNodes;
kNN_struct.NodeTable = NodeTable;
kNN_struct.EdgeTable = EdgeTable;



