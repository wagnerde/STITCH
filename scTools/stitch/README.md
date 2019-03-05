STITCH
=========

STITCH assembles a graph manifold from time series single-cell RNA sequencing data.  STITCH acknowledges that the gene expression signatures defining a complex biological process (e.g. embryonic development) can change dramatically over time.  Rather than projecting data from several timepoints into a single low-dimensional latent space, STITCH instead assembles a manifold that is defined through a series of overlapping, locally-defined PCA subspaces.

1. Non-mutual k-nearest-neighborhoods are first obtained for each cell in timepoint _i_. Neighbor edges are queried from timepoints _i_ (within-timepoint edges) and _i-1_ (link edges) after projecting into a PCA subspace defined by _all_ cells from timepoint _i_.  
2. Outgoing edges are then subject to local and global neighborhood restictions.
3. The graph is restricted to mutual edges.

STITCH is heuristic but it does not use pre-defined clusters to guide construction of the manifold.  Nor does it make topological assumptions (e.g. a "tree-like" organization) about the underlying data.  A STITCH graph can be used as a starting point for data visualization and exploration by generating a force-directed layout (e.g. ForceAtlas2 in [Gephi](https://gephi.org/)).

This version of STITCH was written in Matlab (2017a) and requires the 'Statistics and Machine Learning Toolbox'.

STITCH was developed and used to analyze zebrafish embryonic development in:  
**Single-cell mapping of gene expression landscapes and lineage in the zebrafish embryo.**  Wagner DE, Weinreb C, Collins ZM, Briggs JA, Megason SG, Klein AM. Science 26 Apr 2018. [doi:10.1126/science.aar4362](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362)


## Usage ##

### Inputs ###
Input data are supplied to STITCH as a Matlab object called "DataSet".  "DataSet" is a structure array with the following fields and one record for each timepoint sample:
1. **DataSet.name**. A unique string identifier for each sample timepoint (REQUIRED).
2. **DataSet.ind**.  A numeric index indicating timeseries order (REQUIRED).
3. **DataSet.X**.  A UMI-filtered counts matrix; rows=transcripts, columns=cells (REQUIRED).  X can be supplied as a sparse matrix.  
4. **DataSet.gene_ind**.  An array of gene row indices, e.g. highly variable genes (OPTIONAL). If not provided, variable genes can be determined  using a corrected Fano factor test statistic as reported in [Klein et. al. 2015](https://doi.org/10.1016/j.cell.2015.04.044).
5. **DataSet.batch_flag**. A numeric array of sample batch IDs for each cell within each timepoint (OPTIONAL). If present, gene normalizations are performed within each batch.
6. **DataSet.nDim**. Number of PCA dimensions to use for each timepoint (OPTIONAL). If not provided, a single value ('nDim') will be used for all timepoints.

An example DataSet from [Wagner et. al. 2018](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362) can be downloaded [here](https://kleintools.hms.harvard.edu/paper_websites/wagner_zebrafish_timecourse2018/WagnerScience2018.mat).


### Settings ###
Parameter settings for data preprocessing (see "scTools") and the main STITCH pipeline (see "stitch_get_graph") are specified using optional name/value pairs.  Default behavior implements settings used in [Wagner et. al. 2018](http://science.sciencemag.org/content/early/2018/04/25/science.aar4362).

```
 'k_initial'
           Initial number of nearest neighbors to consider when
           constructing kNN graphs and timepoint links (default=200).

 'k_final'
           Number of nearest neighbors to retain for each node in the
           the final graph (default=20).

 'nDim'
           Number of principal component dimensions to use for subspace
           projections (default=200).

 'distance_metric'
           Distance metric supplied to knnsearch (default='correlation') 

 'graph_max_D_local'
           Local distance threshold for filtering graph edges. Expressed 
           as a multiple of the nearest neighbor distance from the source 
           node (default=3).

 'link_max_D_local'
           Local distance threshold for filtering link edges. Expressed 
           as a multiple of the nearest neighbor distance from the source 
           node (default=3).

 'graph_max_D_global':
           Global distance threshold for filtering graph edges. Expressed 
           as a z-score over all edges in the current timepoint 
           (default=-1).
           
 'link_max_D_global':
           Global distance threshold for filtering link edges. Expressed 
           as a z-score over all edges in the current timepoint 
           (default=0).           
           
 'require_mutual_edges'
           Only retain edges for which both nodes appear in the other 
           node's outgoing edge list (default=true).
```


### Run STITCH ###

We provide a script: "stitch_Wagner2018.m" for running the entire pipeline.  The script begins by downloading data, performing total counts normalization, identifying variable genes, calculating and plotting a STITCH graph.

Run the script by typing the following into the Matlab command line:     
  ```
  run('script_runSTITCH_Wagner2018.m')
  ```

### Visualizing STITCH Graphs ###

We recommend visualizing STITCH graphs interactively using the ForceAtlas2 layout in [Gephi](https://gephi.org/). The Matlab graph object can be imported into Gephi from DOT format, for example by using the 'graph_to_dot' function from [Matlab-Graphviz interface](https://www.mathworks.com/matlabcentral/fileexchange/4518-matlab-graphviz-interface) by Leon Peshkin.  
  ```
  addpath('scTools/graphViz2Mat1')
  graph_to_dot(adjacency(G), 'directed', 0, 'filename', 'Wagner2018_stitch.dot')
  ```

Gephi coordinates can be exported as a .NET file and imported back to Matlab.  The example dataset provides pre-computed XY coordinates.  
  ```
  XY = import_gephi_xy('Wagner2018_stitch.net')
  ```

Plot the STITCH graph using Gephi coordinates. Default behavior colors nodes by timepoint.    
  ```
  stitch_plot_graph(G, XY)
  ```

Alternatively, color nodes by a specific gene.  
```
  stitch_plot_graph(G, XY, DataSet, gene_names_all, 'nodes', 'nanog')
  stitch_plot_graph(G, XY, DataSet, gene_names_all, 'nodes', 'epcam')
  stitch_plot_graph(G, XY, DataSet, gene_names_all, 'nodes', 'msgn1')
  stitch_plot_graph(G, XY, DataSet, gene_names_all, 'nodes', 'sox19a')
```
See below for a full list of stitch_plot_graph input options.
```
INPUTS:
 G               A STITCH graph object
 XY              Matrix of XY coordinates for each node in G 
 DataSet         A STITCH dataset object 
                 (only required if plotting gene expression values)
 gene_names_all  Gene names: a cell array of strings 
                 (only required if plotting gene expression values)

OPTIONAL INPUTS:
'nodes'
                Input string, specifying node format style.  
                 'timepoint': colors nodes by sample/timepoint (default).
                 'none': nodes are hidden from view.
                 'degree'|'closeness'|'eigenvector'|
                 'betweenness'|'pagerank': nodes colored based on graph 
                  centrality scores.
                 Any other input will be interpreted as a gene name and
                  searched against gene_names_all.                 

'node_rgb'       An array of RGB triplet values for coloring all nodes of 
                 the graph (e.g. [255 255 255]).  If specified, this 
                 option overrides the 'nodes' option.

'node_scores' 
                 An array of custom node values, which must match the 
                 number of nodes in the graph.  If specified, this option 
                 overrides the 'nodes' option.

'node_color_scale'
                 Input string: 'log' (default) or 'linear'.

'node_size'
                 Float specifying node marker size (default = 0.8) 

'hide_zero_count_nodes'
                 When plotting gene expression or custom values, option 
                 to hide all nodes with zero counts (default = true).

'edges'
                 Input string, specifying edge format style. 
                 'alpha': thin, semi-transparent edges (default).
                 'black'|'none': edges colored black or hidden from view.
                 'internal'|'bridges': highlight edges that connect nodes
                  within timepoints, or bridging between timepoints.
                 'weights': color edges based on weight (e.g. correlation
                  distance in PCA space).

'edge_scores'
                 An array of custom edge values, which must match the 
                 number of edges in the graph.  If specified, this option 
                 overrides the 'edges' option.

'edge_color_scale'
                 Input string: 'log' or 'linear' (default).

```

### Coarse-Graining the STITCH Graph ###

The STITCH graph may be simplified via a process called 'coarse-graining'.  In the coarse-grained graph, nodes correspond to  clusters of cells (specified by node_IDs), and edges reflect the connectivity of two clusters (the proportion of original shared single-cell edges to total outgoing edges).  Through the coarse-grained graph, a spanning tree "scaffold" is then traced, which highlights the strongest set of edges that bridge clusters between timepoints.

```
  node_IDs = [];
  for j = nTimePoints:-1:1
      node_IDs = [node_IDs; string(DataSet(j).celldata.cell_IDs_names)]; 
  end
  [G_cg, G_cg_scaff] = stitch_coarse_grain(G, index_to_graph(G, node_IDs));
  graph_to_dot(adjacency(G_cg_scaff), 'directed', 0, 'filename', 'gephi/Wagner2018_cg.dot')
```
Plot the Coarse-grained graph, again using imported Gephi coordinates.
```
  XY_cg = import_gephi_xy('gephi/Wagner2018_cg.net');
  figure; stitch_plot_graph_cg(G_cg, XY_cg, [], [], 'node_size', 5)
```

### Exporting STITCH Graphs ###

Export attributes of the STITCH graph as text files (e.g. for import into non-Matlab environments).   
  ```
  stitch_export('export_directory', DataSet, G, XY, gene_names)
  ```

The following files will be exported:
```
genes.txt        Gene names
counts.csv       Genes x cells counts matrix
annot.txt        Annotation flags for each cell (DataSet.celldata)
timepoints.txt   Timepoint/sample flags for each cell (DataSet.name)
edges.csv        Edge table for G
coordinates.txt  XY coordinates for each node of G
```