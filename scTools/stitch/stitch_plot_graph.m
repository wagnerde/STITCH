function stitch_plot_graph(G, XY, DataSet, gene_names_all, varargin)
% Usage: stitch_plot_graph(G, XY, DataSet, gene_names_all, varargin)
%
% Plots a 2D layout of a STITCH graph.  Nodes and edges are formatted
% according to optional name/value pairs.  Default behavior colors nodes by
% sample/timepoint.
% 
% INPUTS:
% G               A STITCH graph object
% XY              Matrix of XY coordinates for each node in G 
% DataSet         A STITCH dataset object (only required if plotting gene
%                 expression values).
% gene_names_all  Gene names: a cell array of strings (only required if 
%                 plotting gene expression values). 
%
% Optional input name/value pairs:
% 'nodes'
%                 Input string, specifying node format style.  
%                  'timepoint': colors nodes by sample/timepoint (default).
%                  'none': nodes are hidden from view.
%                  'degree'|'closeness'|'eigenvector'|
%                  'betweenness'|'pagerank': nodes colored based on graph 
%                   centrality scores.
%                  Any other input will be interpreted as a gene name and
%                   searched against gene_names_all.                 
%
% 'node_rgb'       An array of RGB triplet values for coloring all nodes of 
%                  the graph (e.g. [255 255 255]).  If specified, this 
%                  option overrides the 'nodes' option.
%
% 'node_scores' 
%                  An array of custom node values, which must match the 
%                  number of nodes in the graph.  If specified, this option 
%                  overrides the 'nodes' option.
%
% 'node_color_scale'
%                  Input string: 'log' (default) or 'linear'.
%
% 'node_size'
%                  Float specifying node marker size (default = 0.8) 
% 
% 'hide_zero_count_nodes'
%                  When plotting gene expression or custom values, option 
%                  to hide all nodes with zero counts (default = true).
%
% 'edges'
%                  Input string, specifying edge format style. 
%                  'alpha': thin, semi-transparent edges (default).
%                  'black'|'none': edges colored black or hidden from view.
%                  'internal'|'bridges': highlight edges that connect nodes
%                   within timepoints, or bridging between timepoints.
%                  'weights': color edges based on weight (e.g. correlation
%                   distance in PCA space).
%
% 'edge_scores'
%                  An array of custom edge values, which must match the 
%                  number of edges in the graph.  If specified, this option 
%                  overrides the 'edges' option.
%
% 'edge_color_scale'
%                  Input string: 'log' or 'linear' (default).
%

%% PARAMETER SETTINGS
% Set defaults
def.nodes = 'timepoint';
def.node_rgb = [];
def.node_color_scale = 'log';
def.node_size = 0.8;
def.hide_zero_count_nodes = true;
def.edges = 'alpha';
def.edge_color_scale = 'linear';
def.node_scores = [];
def.edge_scores = []; 

% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'stitch_plot_graph';
parserObj.StructExpand = false; 
parserObj.addOptional('nodes',def.nodes);
parserObj.addOptional('node_rgb',def.node_rgb);
parserObj.addOptional('node_color_scale',def.node_color_scale);
parserObj.addOptional('node_size',def.node_size);
parserObj.addOptional('hide_zero_count_nodes',def.hide_zero_count_nodes);
parserObj.addOptional('edges',def.edges);
parserObj.addOptional('edge_color_scale',def.edge_color_scale);
parserObj.addOptional('node_scores',def.node_scores);
parserObj.addOptional('edge_scores',def.edge_scores);

% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;


%% CODE
% make the base figure
% figure
set(gca,'Position',[0.05 0.05 0.9 0.9])
set(gcf,'Color','[0.98 0.98 0.98]');
p = plot(G, 'XData', XY(:,1), 'YData', XY(:,2));
p.MarkerSize = settings.node_size;
axis off

% create empty containers for node and edge data
NodeCData_tmp = [];
EdgeCData_tmp = [];
    
% custom node/edge scores or rgb values override other options
if ~isempty(settings.node_scores)
    settings.nodes = 'scores';
end
if ~isempty(settings.node_rgb)
    settings.nodes = 'rgb';
end
if ~isempty(settings.edge_scores)
    settings.edges = 'scores';
end

% format nodes
switch settings.nodes
    
    case 'timepoint' % COLOR NODES BY TIMEPOINT
        NodeCData_tmp = grp2idx(G.Nodes.OriginalDataSet);
        colormap jet
        settings.node_color_scale = 'linear';
        title('Timepoints')       
        
    case 'none' % HIDE NODES
        set(p, 'NodeColor', 'none')
        
    case 'degree'
        p.NodeColor = 'black';
        NodeCData_tmp = degree(G);
        colormap jet
        title('Centrality: Degree')
        
    case 'closeness'
        p.NodeColor = 'black';
        NodeCData_tmp = centrality(G,'closeness');
        colormap jet
        title('Centrality: Closeness')
            
    case 'eigenvector'
        p.NodeColor = 'black';
        NodeCData_tmp = centrality(G,'eigenvector');
        colormap jet
        title('Centrality: Eigenvector')
     
    case 'betweenness'
        p.NodeColor = 'black';
        node_betweenness = centrality(G,'betweenness');
        NodeCData_tmp = node_betweenness;
        colormap jet
        title('Centrality: Betweenness')
        
    case 'pagerank'
        p.NodeColor = 'black';
        NodeCData_tmp = centrality(G,'pagerank');
        colormap jet
        title('Centrality: Pagerank')
    
    case 'rgb'  
        if any(settings.node_rgb>1)
            settings.node_rgb = settings.node_rgb / 255;
        end
        p.NodeColor = settings.node_rgb; 
        
    case 'scores' % COLOR NODES ACCORDING TO CUSTOM SCORES
        NodeCData_tmp = settings.node_scores;      
        if settings.hide_zero_count_nodes
            nodes_with_zero_counts = find(NodeCData_tmp == 0);
            highlight(p, nodes_with_zero_counts, 'MarkerSize', 0.01)  
        end
       
    otherwise % COLOR NODES BY COUNTS OF A SPECIFIC GENE
        gene_ind = strcmp(gene_names_all, settings.nodes);
        cell_ind = str2num(cell2mat(G.Nodes.Name));
        NodeCData_tmp = [];
        X_array = {DataSet.X_norm};
        for j=length(X_array):-1:1 % concatenate rows for this gene on the fly...
            NodeCData_tmp = [NodeCData_tmp full(X_array{j}(gene_ind,:))];
        end        
        NodeCData_tmp = NodeCData_tmp(cell_ind);
        colormap jet
        if settings.hide_zero_count_nodes 
            nodes_with_zero_counts = find(NodeCData_tmp == 0);
            highlight(p, nodes_with_zero_counts, 'MarkerSize', 0.01)  
        end
        title(['Counts: ' settings.nodes])
        
end


% format edges
switch settings.edges
    
    case 'black' % BLACK EDGES
        set(p, 'EdgeColor', 'black', 'EdgeAlpha', 0.5)
    
    case 'alpha' % BLACK EDGES WITH TRANSPARENCY
        set(p, 'EdgeColor', 'black', 'EdgeAlpha', 0.05)

    case 'none' % HIDE EDGES    
        set(p, 'EdgeColor', 'none')
        
    case 'internal' % HIGHLIGHT INTERNAL TIMEPOINT EDGES
        EdgeCData_tmp = G.Edges.InternalEdge;
        colormap([[1 1 1];[0 0 0]])
        title('Internal Timepoint Edges')
    
    case 'bridges' % HIGHLIGHT EDGES CONNECTING TIMEPOINTS
        EdgeCData_tmp = G.Edges.LinkEdge;
        colormap([[1 1 1];[0 0 0]])
        title('Internal Timepoint Edges')
        title('Timepoint Bridges')
       
    case 'weights' % COLOR EDGES ACCORDING TO THEIR WEIGHTS
        EdgeCData_tmp = G.Edges.D_orig;
        colormap(jet)
        title('Edge Weights')
        
    case 'scores' % COLOR EDGES ACCORDING TO CUSTOM SCORES
        set(p, 'EdgeColor', 'black', 'EdgeAlpha', 0.2, 'Linewidth', 0.5)
        EdgeCData_tmp = settings.edge_custom;
        colormap((jet))

end

% apply formats to plot
if ~isempty(NodeCData_tmp)
    switch settings.node_color_scale
        case 'log'
            p.NodeCData = log10(NodeCData_tmp);
        case 'linear'
            p.NodeCData = NodeCData_tmp;         
    end
end

if ~isempty(EdgeCData_tmp)
    switch settings.edge_color_scale
        case 'log'
            p.EdgeCData = log10(EdgeCData_tmp);
        case 'linear'
            p.EdgeCData = EdgeCData_tmp;
    end
end

caxis('auto')

