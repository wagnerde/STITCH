function kNN_struct = get_knn(X_norm, gene_ind, nDim, k_use, distance_metric, coeff)
% Usage: kNN_struct = get_knn(X_norm, gene_ind, nDim, k_use, distance_metric, coeff)
%
% Identify k-nearest neighbors from a normalized single-cell counts matrix.
% Begins by reducing dimensionality w/ PCA, identifying neighbors with
% knnsearch, and then returning node/edge tables. If supplied, precomputed 
% PCA coefficients can also be used.
%
% INPUTS:
% 'X_norm'      
%         A total-counts normalized single-cell UMI counts matrix
%         rows=transcripts, columns=cells
%
% 'gene_ind'    
%         Row indices for genes of interest (e.g. highly variable genes)  
%
% 'nDim'
%         Number of PCA dimensions to use.
%
% 'k_use' 
%         Number of nearest neighbors to retain for each cell.  
%
% 'distance_metric' 
%         Distance metric to supply to knnsearch (e.g. 'euclidean').
% 
% 'coeff'
%         Optional pre-computed PCA coefficients.
%
% 
% OUTPUTS:
% kNN_struct 
%         Structure array containing node tables and edge tables
%
%

%% CODE

% subset and standardize data matrix
nNodes = size(X_norm,2);
X_indexed = X_norm(gene_ind,:);
X_z = zscore(X_indexed')';

% perform PCA
if isempty(coeff)
    [~,scores,~,~,~,~] = pca(X_z', 'Centered', true, 'Economy', true, 'NumComponents', nDim);
    scores(:,nDim+1:end) = [];
else
    scores = project_to_PCs(X_z, coeff, nDim);
    scores(:,nDim+1:end) = [];
end

% search for neighbors
k_use = k_use + 1;
[IDX_node2, D] = knnsearch(scores, scores, 'K', k_use, 'Distance', distance_metric);
IDX_node1 = repmat((1:(nNodes))',1,k_use);

% record the nearest neighbor distances for plotting later (column1=self-edges)
D_shortest = D(:,2);

% re-scale D as a multiple of the nearest neighbor distance
D_shortest_rep = repmat(D_shortest,1,k_use);
D_local = D ./ D_shortest_rep;

% convert matrices to 1d arrays; omit column1 (self-edges)
IDX_node2 = mat2array(IDX_node2(:,2:end));
IDX_node1 = mat2array(IDX_node1(:,2:end));
D = mat2array(D(:,2:end));
D_local = mat2array(D_local(:,2:end));

% normalize edge distance to a global z-score
D_global = zscore(D);

% generate labels   
nEdges = length(IDX_node2);
NodeIDs = (1:nNodes)';

% assemble edge & node tables
EdgeTable = table([IDX_node1, IDX_node2], D_global, D_local, D, 'VariableNames', {'EndNodes', 'Weights', 'D_local', 'D_orig'});
NodeTable = table(NodeIDs, 'VariableNames', {'Name'});
NodeTable.Name = cellstr(num2str(NodeTable.Name));

% outputs
kNN_struct.nEdges = nEdges;
kNN_struct.nNodes = nNodes;
kNN_struct.NodeTable = NodeTable;
kNN_struct.EdgeTable = EdgeTable;



