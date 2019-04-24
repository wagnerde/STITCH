function DataSet = stitch_get_louvain_clusters(DataSet, varargin)
% Usage: DataSet = stitch_get_louvain_clusters(DataSet, varargin)
%
% Identifies Louvain cell clusters for each timepoint entry of DataSet.  
% Returns cluster assignments as a column in the table: DataSet.celldata
%

%% PARAMETER SETTINGS
% Set defaults
def.k_force = [];
def.distance_metric = 'correlation';
def.resolution = 1;

% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'stitch_get_Louvain_clusters';
parserObj.StructExpand = false; 
parserObj.addOptional('k_force',def.k_force);
parserObj.addOptional('distance_metric',def.distance_metric);
parserObj.addOptional('resolution',def.resolution);

% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% CODE:

base_counter = 0;

for j = 1:length(DataSet)
    
    disp(['Clustering timepoint ' num2str(j) ' ...'])
    
    % check to see if k_force was specified, otherwise, k_use = sqrt(n_ cells)
    if isempty(settings.k_force)
        k_use = round(sqrt(size(DataSet(j).X,2)));
    else
        k_use = settings.k_force;
    end
    
    % build a knn graph
    tmp.tables = get_knn(full(DataSet(j).X_norm), DataSet(j).gene_ind, DataSet(j).nDim, k_use, settings.distance_metric, []);
    tmp.tables.EdgeTable = filter_duplicate_edges(tmp.tables.EdgeTable);    
    tmp.G = graph(tmp.tables.EdgeTable, tmp.tables.NodeTable);
    
    % perform gen-louvain clustering
    DataSet(j).celldata = table(get_louvain_clusters(tmp.G, settings.resolution)+base_counter, 'VariableNames', {'Louvain_ID'});
    
    % update cluster ID counter
    base_counter = max(DataSet(j).celldata.Louvain_ID);
    
end