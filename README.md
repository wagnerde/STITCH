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
4. **DataSet.gene_ind**.  An array of gene row indices, e.g. highly variable genes (OPTIONAL). If not provided, variable genes will be determined  using a corrected Fano factor test statistic as reported in [Klein et. al. 2015](https://doi.org/10.1016/j.cell.2015.04.044).
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

We provide a script: "stitch_Wagner2018.m" for running the entire pipeline, starting with raw UMI counts tables.  The script begins by downloading counts data, performing total counts normalization, identifying variable genes, calculating and plotting a STITCH graph.

Run the script by typing the following into the Matlab command line:     
  ```run('stitch_Wagner2018.m')```;

We recommend visualizing STITCH graphs using the ForceAtlas2 layout in [Gephi](https://gephi.org/). The Matlab graph object can be imported into Gephi from DOT format, for example by using the 'graph_to_dot' function from [Matlab-Graphviz interface](https://www.mathworks.com/matlabcentral/fileexchange/4518-matlab-graphviz-interface) by Leon Peshkin.  
  ```graph_to_dot(adjacency(G), 'directed', 0, 'filename', 'gephi_graph.dot')```

Gephi coordinates can be exported as a .NET file and imported back to Matlab.  The example dataset provides pre-computed XY coordinates.  
  ```XY = stitch_import_gephi('gephi_export.net')```

Plot the STITCH graph using Gephi coordinates.    
  ```stitch_plot_graph(G, XY)```












