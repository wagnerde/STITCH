function gephi_XY = import_gephi_xy(filename)
% Usage: gephi_XY = import_gephi_xy(filename)
%
% Imports x-y coordinates from Gephi. 
% Input file is a Gephi graph exported in .NET format
% gephi_XY is a matrix of plotting coordinates for each node.
%
%% CODE:
% initialize textscan variables
delimiter = ' ';
startRow = 2;
formatSpec = '%q%q%q%q%q';

% open file
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

% find the row where node data ends
ind_NodeDataEnd = find(strcmp('*Edges',dataArray{1}))-1;
if isempty(ind_NodeDataEnd)
    try
        ind_NodeDataEnd = find(strcmp('*Arcs',dataArray{1}))-1;
    end
end

% read in names, x-coordinates, and y-coordinates for each node
gephi_Id = cellfun(@str2num, dataArray{2}(1:ind_NodeDataEnd));
gephi_XY(:,1) = cellfun(@str2num, dataArray{3}(1:ind_NodeDataEnd));
gephi_XY(:,2) = cellfun(@str2num, dataArray{4}(1:ind_NodeDataEnd));
[~,ind] = sort(gephi_Id);

% export
gephi_XY = gephi_XY(ind,:);


