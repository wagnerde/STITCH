function plot_knn_graph(G, XY, NodeCData)
% Usage: plot_knn_graph(G, XY, NodeCData)
%
% Plots a 2D layout of a knn graph.
% 
% INPUTS:
% G          A Matlab graph object  
% XY         Matrix of XY coordinates for each node in G 
% NodeCData  
%

%% CODE

figure('Position', [500 500 500 500])
set(gca,'Position',[0.05 0.05 0.9 0.9])

p = plot(G, 'XData', XY(:,1),'YData', XY(:,2));
axis square; axis off

p.MarkerSize = 3;
p.EdgeAlpha = 0.2;
p.EdgeColor = [0 0 0];
p.LineWidth = 0.05;
p.NodeColor = [0 0 0];

if exist('NodeCData', 'var')
    p.NodeCData = log10(NodeCData);
    colormap(jet)
end


