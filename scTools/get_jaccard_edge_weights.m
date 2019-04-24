function [A_jw, jaccard_similarities] = get_jaccard_edge_weights(G)
% Usage:
%
%
%
%% CODE:

% get graph nodes (unique names)
unique_node_names = str2num((cell2mat(G.Nodes.Name)));
nNodes = length(unique_node_names);
% get graph edges (unique names)
unique_endnode_names1 = str2num(cell2mat(G.Edges.EndNodes(:,1)));
unique_endnode_names2 = str2num(cell2mat(G.Edges.EndNodes(:,2)));
nEdges = length(unique_endnode_names1);
% get current indices for each endnode
[~,EndNode1] = ismember(unique_endnode_names1,unique_node_names);
[~,EndNode2] = ismember(unique_endnode_names2,unique_node_names);

% get neighborhood for each node
for j = 1:nNodes
    neighborhoods{j} = neighbors(G, j);
end

% get jaccard similarity between node neighborhoods for each edge
for j = 1:nEdges
    jaccard_similarities(j) = length(intersect(neighborhoods{EndNode1(j)},neighborhoods{EndNode2(j)})) / length(union(neighborhoods{EndNode1(j)},neighborhoods{EndNode2(j)}));
end

% generate weighted adjacency matrix
[s,t] = findedge(G);
A_jw = sparse(s, t, jaccard_similarities, nNodes, nNodes);
A_jw = A_jw + A_jw';


