function [array, ind] = index_to_graph(G, array)
%
% Indexes an array to the node names contained in graph G
%
%% CODE:

ind = str2num(cell2mat(G.Nodes.Name));
array = array(ind);


