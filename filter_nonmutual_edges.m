function EdgeTable = filter_nonmutual_edges(EdgeTable)
% Usage: EdgeTable = filter_nonmutual_edges(EdgeTable)
%
% Accepts a table of STITCH edges. 
% Identifies non-duplicate (non-mutual) edges and discards them.
%

%% CODE:

% sort edges so that the higher node index appears in the 2nd of two columns
[~,ia,ic] = unique(sort(EdgeTable.EndNodes(:,[1,2]), 2 ), 'rows');    

% mutual edges should now be duplicate rows, non-mutual edges will have only a single row
nOccurences = histc(ic,unique(ic));

% remove non-mutual edges
is_nonmutual_edge = (nOccurences == 1); 
ind_nonmutual_edge = ia(is_nonmutual_edge);
EdgeTable(ind_nonmutual_edge,:) = [];
