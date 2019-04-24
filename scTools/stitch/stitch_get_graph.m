function G = stitch_get_graph(DataSet, varargin)
% Usage: G = stitch_get_graph(DataSet, varargin)
%
% Main function for the STITCH pipeline.  Inputs a DataSet object and 
% performs the following steps:
% (1) Constructs a kNN graph for each individual timepoint
% (2) Constructs a table of single-cell edges linking adjacent timepoints
% (3) Merges node and edge tables from all timepoints
% (4) Filters to remove weak and non-mutual edges 
%
% INPUTS:
% DataSet  
%           Structure array with multiple fields, and one record for each
%           timepoint. Timepoints should be ordered sequentially.
%           DataSet.ind         is a numeric index of the current timepoint. (required)
%           DataSet.name        is a unique string identifier. (required) 
%           DataSet.X           is a single-cell counts matrix; (required)   
%                               rows=transcripts, columns=cells.  
%           DataSet.gene_ind    is an array of gene row indices, e.g. 
%                               highly-variable genes. (optional)
%           DataSet.batch_flag  is a numeric index of individual sample 
%                               batches within a timepoint. if present, 
%                               gene normalizations are performed within 
%                               each batch. (optional)
%           DataSet.nDim        indicates the number of PCA dimensions to 
%                               to use for each timepoint (optional).
%                               If not provided, a single value ('nDIM')
%                               will be used for all timepoints.
%                               
%                               
% Optional input name/value pairs: 
% 'k_initial'
%           Initial number of nearest neighbors to consider when
%           constructing kNN graphs and timepoint links (default=200).
%
% 'k_final'
%           Number of nearest neighbors to retain for each node in the
%           final graph (default=20).
%
% 'nDim'
%           Number of principal component dimensions to use for subspace
%           projections (default=200).
%
% 'distance_metric'
%           Distance metric supplied to knnsearch (default='correlation') 
%
% 'graph_max_D_local'
%           Local distance threshold for filtering graph edges. Expressed 
%           as a multiple of the nearest neighbor distance from the source 
%           node (default=3).
%
% 'link_max_D_local'
%           Local distance threshold for filtering link edges. Expressed 
%           as a multiple of the nearest neighbor distance from the source 
%           node (default=3).
%
% 'graph_max_D_global':
%           Global distance threshold for filtering graph edges. Expressed 
%           as a z-score over all edges in the current timepoint 
%           (default=-1).
%           
% 'link_max_D_global':
%           Global distance threshold for filtering link edges. Expressed 
%           as a z-score over all edges in the current timepoint 
%           (default=0).
%
% 'require_mutual_edges'
%           Only retain edges for which both nodes appear in the other 
%           node's outgoing edge list (default=true).
%
% OUTPUTS:
% G         Matlab graph object containing node and edge tables
% DataSet   DataSet object with additional data preprocessing fields added
%           

%% PARAMETER SETTINGS
% Set defaults
def.k_initial = 200;
def.k_final = 20;
def.nDim = 200;
def.distance_metric = 'correlation';
def.graph_max_D_local = 3;  
def.link_max_D_local =  3; 
def.graph_max_D_global = -1;
def.link_max_D_global = 0; 
def.require_mutual_edges = true;

% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'stitch';
parserObj.StructExpand = false; 
parserObj.addOptional('k_initial',def.k_initial);
parserObj.addOptional('nDim',def.nDim);
parserObj.addOptional('k_final',def.k_final);
parserObj.addOptional('distance_metric',def.distance_metric);
parserObj.addOptional('graph_max_D_local',def.graph_max_D_local);
parserObj.addOptional('link_max_D_local',def.link_max_D_local);
parserObj.addOptional('graph_max_D_global',def.graph_max_D_global);
parserObj.addOptional('link_max_D_global',def.link_max_D_global);
parserObj.addOptional('require_mutual_edges',def.require_mutual_edges);

% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% CODE:

% Check whether preprocessing is still needed
if (~isfield(DataSet, 'X_norm') || ~isfield(DataSet, 'tot_counts') || ~isfield(DataSet, 'gene_ind'))
    disp('Missing Fields.  Please run "get_normalized_counts" and "get_variable_genes" before calculating a STITCH graph.')
    return
end

% Step1. Build kNN graphs within timepoints
disp('Building kNN graphs...')
nTimePoints = length(DataSet);
for j = nTimePoints:-1:1
    stitch_obj.knn(nTimePoints-j+1) = stitch_get_knn(DataSet(j), settings);
end

% Step2. Get edges linking adjacent timepoints
for j = nTimePoints-1:-1:1 
    stitch_obj.links(j) = stitch_get_links(DataSet(j+1), DataSet(j), settings);
end

% Step3. Merge all node and edge tables
stitch_obj.merged = stitch_merge(stitch_obj.knn, stitch_obj.links);

% Step4. Filter non-mutual & weak edges
stitch_obj.filtered = stitch_filter(stitch_obj.merged, settings);

% Export graph
G = stitch_obj.filtered.G;

