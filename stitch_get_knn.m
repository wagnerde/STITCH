function stitch_knn = stitch_get_knn(DataSet, settings)
% Usage: stitch_knn = stitch_get_knn(DataSet, settings)
%
% Generates a kNN graph from single record the DataSet structure array.
% Begin by reducing dimensionality with PCA, identifying nearest neighbors 
% with knnsearch, and returning node and edge tables.
% 
% INPUTS:
% DataSet       A single record/row from the DataSet structure array
% settings      Structure containing STITCH settings
%
% OUTPUT:
% stitch_knn    A structure array containing node and edge tables for a
%               single timepoint
% 

%% CODE:

% if provided, import timepoint-specific nDim setting
if isfield(DataSet,'nDim')
    settings.nDim = DataSet.nDim;
end

% subset data by highly variable genes
nGenes = length(DataSet.gene_ind);
X_indexed = DataSet.X_norm(DataSet.gene_ind,:);
nNodes = size(DataSet.X_norm,2);

% convert transcript counts to zscores
X_indexed = full(X_indexed);
if isfield(DataSet,'batch_flag')
    X_z = get_zscore_batch(X_indexed, DataSet.batch_flag);
else
    X_z = zscore(X_indexed,0,2);
end

% perform PCA
[~,scores] = pca(X_z', 'Centered', true, 'Economy', true, 'NumComponents', min(nGenes, settings.nDim));

% search for neighbors
settings.k_initial = settings.k_initial + 1;
[IDX_node2, D] = knnsearch(scores, scores, 'K', settings.k_initial, 'Distance', settings.distance_metric);
IDX_node1 = repmat((1:(nNodes))',1,settings.k_initial);

% record the nearest neighbor distances for plotting later (column1=self-edges, column2=first neighbor)
D_shortest = D(:,2);

% re-scale D as a multiple of the nearest neighbor distance
D_shortest_rep = repmat(D_shortest,1,settings.k_initial);
D_local = D ./ D_shortest_rep;

% convert matrices to 1d arrays by row; omit column1 (self-edges)
IDX_node2 = mat2array(IDX_node2(:,2:end));
IDX_node1 = mat2array(IDX_node1(:,2:end));
D = mat2array(D(:,2:end));
D_local = mat2array(D_local(:,2:end));

% normalize edge distances by z-score
D_global = zscore(D);

% generate labels   
nEdges = length(IDX_node2);
EdgeLabels = repmat({DataSet.name}, nEdges, 1);
NodeIDs = (1:nNodes)';
NodeLabels = cellstr(repmat(DataSet.name, nNodes, 1));
LinkEdge = zeros(nEdges,1);
InternalEdge = ones(nEdges,1);
OriginalDataSet = repmat(DataSet.ind,nNodes,1);
BasisSpace = repmat(DataSet.ind,nEdges,1);
BasisOutgoing = ones(nNodes*(settings.k_initial-1),1);

% assemble edge & node tables
EdgeTable = table([IDX_node1, IDX_node2], D_global, EdgeLabels, BasisSpace, LinkEdge, InternalEdge, BasisOutgoing, D_global, D_local, D, 'VariableNames', {'EndNodes', 'Weights', 'EdgeLabels', 'BasisSpace', 'LinkEdge', 'InternalEdge', 'BasisOutgoing', 'D_global', 'D_local', 'D_orig'});
NodeTable = table(NodeIDs, NodeLabels, NodeIDs, OriginalDataSet, 'VariableNames', {'Name', 'NodeLabels', 'OriginalName', 'OriginalDataSet'});

% outputs
stitch_knn.set_name = DataSet.name;
stitch_knn.DataSet_ind = DataSet.ind;
stitch_knn.nEdges = nEdges;
stitch_knn.nNodes = nNodes;
stitch_knn.NodeTable = NodeTable;
stitch_knn.EdgeTable = EdgeTable;
stitch_knn.scores = scores;

disp(['Done building kNN graph: ' stitch_knn.set_name]);

%% HELPER FUNCTIONS:

function array = mat2array(mat1)
%
% Converts a matrix to an array, row by row
%
mat1 = mat1';
array = mat1(:);
end

end
