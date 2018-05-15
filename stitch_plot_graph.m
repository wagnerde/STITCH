function stitch_plot_graph(G, XY)
% Usage: stitch_plot_graph(G, XY)
%
% Plots a 2D layout of a STITCH graph.  Nodes are colored by timepoint.
% 
% INPUTS:
% G         A STITCH graph object  
% XY        Matrix of XY coordinates for each node in G 
%

%% CODE


figure
set(gca,'Position',[0.05 0.05 0.9 0.9])
p = plot(G, 'XData',XY(:,1),'YData',XY(:,2), 'NodeCData', grp2idx(G.Nodes.OriginalDataSet));
p.MarkerSize=1;
p.EdgeAlpha = 0.05;
p.EdgeColor = [0 0 0];
colormap(jet)
axis off


