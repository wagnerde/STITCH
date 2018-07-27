function EdgeTable = filter_duplicate_edges(EdgeTable)
% Usage: EdgeTable = filter_duplicate_edges(EdgeTable)
%
% Accepts a table of STITCH edges. 
% Identifies and removes any duplicate edges.
% Non-duplicate edges are retained.
%

%% CODE:

% sort edges so that the higher node index appears in the 2nd of two columns
[~,ia,ic] = unique(sort(EdgeTable.EndNodes(:,[1,2]), 2 ), 'rows');    

% mutual edges should now be duplicate rows, non-mutual edges will have only a single row
nOccurences = histc(ic,unique(ic));

% keep only the first edge from each duplicate pair
is_duplicate_edge = (nOccurences == 2); 
ind_first_edge_from_each_duplicate_pair = ia(is_duplicate_edge);  

% but also keep any non-duplicate edges
is_non_duplicate_edge = (nOccurences == 1); 
ind_non_duplicate_edge = ia(is_non_duplicate_edge);  

% update the edge table
keep_ind = [ind_non_duplicate_edge; ind_first_edge_from_each_duplicate_pair];
EdgeTable = EdgeTable(keep_ind,:);
