function [G_giant, G_others] = filter_graph_giant_component(G)
% Usage: [G_giant, G_others] = filter_graph_giant_component(G)
%
% Identifies and returns the largest connected component in a graph. 
% Inputs and Outputs are all Matlab graph objects.
%
%
%% CODE: 

% identify connected components and the # of nodes in each
node_component_ids = conncomp(G);
[~,~,ic] = unique(node_component_ids);
component_sizes = histc(ic,unique(ic));
component_sizes_top = sort(component_sizes, 'descend');

% find the largest connected component
[~,giant_component_id] = max(component_sizes);
giant_component_members = find(node_component_ids == giant_component_id);
G_giant = subgraph(G, giant_component_members);

% all remaining components
small_component_members = find(node_component_ids ~= giant_component_id);
G_others = subgraph(G, small_component_members);


