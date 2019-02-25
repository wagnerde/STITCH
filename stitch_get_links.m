function stitch_links = stitch_get_links(DataSet_basis, DataSet_link, settings)
% Usage: stitch_links = stitch_get_links(DataSet_basis, DataSet_link, settings)
%
% Generates an edge list linking between timepoints. Each timepoint as 
% inputted as a single row/record from the DataSet structure array.
% Begin by reducing dimensionality w/ PCA using a specified 'basis' matrix. 
% A second matrix is then projected into the PCA subspace of the 'basis' 
% matrix.  The two projected datasets are then used to perform knnsearch, 
% returning a set of node/edge tables. 
%
% INPUTS:
% DataSet_basis   A single record/row from the DataSet structure array,
%                 which is also used to define the PCA subspace for both
%                 timepoints.  Typically this is the later timepoint.
% DataSet_link    A single record/row from the DataSet structure array,
%                 which is projected into the PCA subspace of DataSet_basis.  
%                 Typically this is the earlier timepoint.
% settings        Structure containing STITCH settings
%
% OUTPUT:
% stitch_knn    A structure array containing node and edge tables for a
%               timepoint-timepoint link.

%% CODE:
if DataSet_basis.ind == DataSet_link.ind
    disp('Error! To map within a single timepoint, use "stitch_get_knn" instead')
    return
end

% if provided, import timepoint-specific nDim setting
if isfield(DataSet_basis,'nDim')
    settings.nDim = DataSet_basis.nDim;
end

% determine whether the basis dataset occurs earlier or later in time
if DataSet_basis.ind < DataSet_link.ind 
    basis_first = true;
else
    basis_first = false;
end

% subset both data matrices by basis gene indices
nGenes = length(DataSet_basis.gene_ind);
X_basis_indexed = DataSet_basis.X_norm(DataSet_basis.gene_ind,:);
X_link_indexed = DataSet_link.X_norm(DataSet_basis.gene_ind,:);

% convert transcript counts to zscores
X_basis_indexed = full(X_basis_indexed);
X_link_indexed = full(X_link_indexed);
if isfield(DataSet_basis,'batch_flag')
    X_basis_z = get_zscore_batch(X_basis_indexed, DataSet_basis.batch_flag);
    X_link_z = get_zscore_batch(X_link_indexed, DataSet_link.batch_flag);
else
    X_basis_z = zscore(X_basis_indexed, 0, 2);
    X_link_z = zscore(X_link_indexed, 0, 2);
end

% perform PCA on the basis dataset
[coeff_basis,scores_basis] = pca(X_basis_z', 'Centered', true, 'Economy', true, 'NumComponents', min(nGenes, settings.nDim));

% project the link dataset into basis dataset PC subspace
scores_projection = X_link_z'*coeff_basis;

% merge the two datasets, basis first
scores_merged = [scores_basis; scores_projection];

% search for neighbors within and between the two datasets
nNodes_basis = size(X_basis_indexed,2);
nNodes_projection = size(X_link_indexed,2);
settings.k_initial = settings.k_initial + 1;
[IDX_node2, D] = knnsearch(scores_merged, scores_merged, 'K', settings.k_initial, 'Distance', settings.distance_metric);
IDX_node1 = repmat((1:(nNodes_basis+nNodes_projection))',1,settings.k_initial);

% retain the nearest neighbor distances (column1=self-edges, column2=first neighbor)
D_shortest = D(:,2);

% re-scale D as a multiple of the nearest neighbor distance ('D_local')
D_shortest_rep = repmat(D_shortest,1,settings.k_initial);
D_local = D ./ D_shortest_rep;

% convert matrices to 1d arrays by row; omit column1 (self-edges)
IDX_node2 = mat2array(IDX_node2(:,2:end));
IDX_node1 = mat2array(IDX_node1(:,2:end));
D = mat2array(D(:,2:end));
D_local = mat2array(D_local(:,2:end));

% flag edges that link the two datasets
% case 1: edges that start in projection and connect to basis
% case 2: edges that start in basis and connect to projection
boundary_ind = nNodes_basis;
case1_ind = (IDX_node1 <= boundary_ind) & (IDX_node2 > boundary_ind);        
case2_ind = (IDX_node1 > boundary_ind) & (IDX_node2 <= boundary_ind);
ind_link_edges = case1_ind | case2_ind;

% normalize edge distance by z-score
D_global = zscore(D);

% generate labels       
set_name = [DataSet_link.name '->' DataSet_basis.name];
nEdges = length(IDX_node1);
EdgeLabels = repmat({set_name}, nEdges, 1);   
LinkEdge = ones(nEdges,1);
InternalEdge = zeros(nEdges,1); 
BasisSpace = repmat(DataSet_basis.ind,nEdges,1);
BasisOutgoing = [ones(nNodes_basis*(settings.k_initial-1),1); zeros(nNodes_projection*(settings.k_initial-1),1)];

% assemble edge table
EdgeTable = table([IDX_node1, IDX_node2], D_global, EdgeLabels, BasisSpace, LinkEdge, InternalEdge, BasisOutgoing, D_global, D_local, D, 'VariableNames', {'EndNodes', 'Weights', 'EdgeLabels', 'BasisSpace', 'LinkEdge', 'InternalEdge', 'BasisOutgoing', 'D_global', 'D_local', 'D_orig'});

% keep only the between-timepoint edges
EdgeTable = EdgeTable(ind_link_edges,:);

% outputs
stitch_links.set_name = set_name;
stitch_links.basis_ind = DataSet_basis.ind;
stitch_links.link_ind = DataSet_link.ind;
stitch_links.basis_first = basis_first;
stitch_links.nEdges = height(EdgeTable);
stitch_links.EdgeTable = EdgeTable;

disp(['Done stitching: ' stitch_links.set_name]);

%% HELPER FUNCTIONS:

function array = mat2array(mat)
%
% Converts a matrix to an array, row by row
%
mat = mat';
array = mat(:);
end




end

