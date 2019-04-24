function spring_export(path, DataSet, G, XY, gene_names)
% Usage: spring_export(path, DataSet, G, XY, gene_names)
%
% Export attributes of a SPRING-like graph (e.g. for import into ScanPy):
% 
% INPUTS:
% path              Name of directory to write files (string).
% DataSet           Single entry of DataSet
% G                 SPRING-like Graph object
% XY                XY coordinates imported from Gephi
% gene_names        Object containing gene names (cell array of strings)
%
% OUTPUT FILES:
%  genes.txt        Gene names
%  counts.csv       Genes x cells counts matrix
%  annot.txt        Annotation flags for each cell (DataSet.celldata)
%  timepoints.txt   Timepoint/sample flags for each cell (DataSet.name)
%  edges.csv        Edge table for G
%  coordinates.txt  XY coordinates for each node of graph object G
%
%% SETTINGS:
%path = 'export';

%% CODE:

% check that export directory exists
if ~exist(path, 'dir')
    mkdir(path);
end

% export gene names
disp('Exporting gene names...')
writetable(cell2table(gene_names), [path '/' 'genes.txt'],'WriteVariableNames',0)

% get original indices of cells in G
tmp_cell_ind = str2num(cell2mat(G.Nodes.Name));

% export counts matrix
% disp('Exporting counts matrix...')
% tmp_X = [];
% nTimePoints = length(DataSet);
% for j = nTimePoints:-1:1
%     tmp_X = [tmp_X DataSet(j).X]; 
% end
% tmp_X = tmp_X(:,tmp_cell_ind);
% dlmwrite([path '/' 'counts.csv'], full(tmp_X)',',')
% clear tmp_X

% export flags for each annotation column in DataSet.celldata
disp('Exporting annotations...')
nAnnots = length(DataSet.celldata.Properties.VariableNames);
for k = 1:nAnnots
    tmp_IDs = DataSet.celldata(:,k); 
    tmp_IDs = tmp_IDs(tmp_cell_ind,:);
    writetable(tmp_IDs, [path '/' DataSet.celldata.Properties.VariableNames{k} '.txt'],'WriteVariableNames',0);
end

% edge table
disp('Exporting edge list...')
dlmwrite([path '/' 'edges.csv'], get_indexed_edge_list(G), ';');

% export XY coordinates
disp('Exporting XY coordinates...')
tmp_ind = 0:length(XY)-1;
tmp_XY = [tmp_ind' XY(:,1) -1*XY(:,2)]; % invert Y axis to match SPRING layout
dlmwrite([path '/' 'coordinates.txt'],tmp_XY,',')
disp('Done!')
