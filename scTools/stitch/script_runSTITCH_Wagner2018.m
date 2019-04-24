%% SCRIPT: Generate STITCH graph from zebrafish timecourse scRNAseq data
%
% This script contains all the necessary steps for performing total counts
% normalization, identifying variable genes, calculating and plotting the 
% Wagner 2018 STITCH and coarse-grained graphs.
% 

%% (1) Download Wagner 2018 counts data
if ~(exist('WagnerScience2018.mat', 'file') == 2)
    disp('Downloading data...')
    websave('WagnerScience2018.mat','https://kleintools.hms.harvard.edu/paper_websites/wagner_zebrafish_timecourse2018/WagnerScience2018.mat');
end

%% (2) Load mat file
disp('Opening mat file...')
load('WagnerScience2018.mat')

%% (3) Add scTools functions to Matlab path
addpath('scTools')

%% (4) Perform total counts normalization 
nTimePoints = length(DataSet);
disp('Normalizing data...')
for j = 1:nTimePoints
    [DataSet(j).X_norm, DataSet(j).tot_counts] = get_normalized_counts(DataSet(j).X);
end
disp('Done normalizing')

%% (5) Identify variable genes
% Define settings used in Wagner et al. 2018.  
disp('Getting variable genes...')
CV_eff =   [0.60 0.40 0.40 0.40 0.40 0.40 0.40];
CV_input = [0.20 0.15 0.15 0.20 0.22 0.22 0.22];
exclude_genes = {'cdk1','mcm2', 'mcm7', 'rrm2', 'cenpa', 'cdc6', 'ccnf', 'cdca4', ...
                 'ccnd1', 'kif4','hmgb1b', 'hmgb3a', 'hspd1', 'hspa9', 'rplp0', ...
                 'hnrnpaba', 'rps2', 'rps12', 'rpl12', 'rps13', 'rps14', 'rps15a', ...
                 'rpl10', 'rps3a', 'rpl31', 'rpl37', 'rps6', 'rpl9', 'rpl11', ...
                 'rpl34', 'rpl13', 'rpl36a', 'rpl26', 'rps8a', 'rpl21', 'rps27.1',...
                 'rpl27a', 'cirbpb'};
for j = 1:nTimePoints
    DataSet(j).gene_ind = get_variable_genes(DataSet(j), 'CV_eff', CV_eff(j) , 'CV_input', CV_input(j), 'excludeGeneNames', exclude_genes, 'allGeneNames', gene_names_all);
end
disp('Done identifying variable genes.')

%% (6) Calculate the STITCH graph
disp('Calculating STITCH graph...')
G = stitch_get_graph(DataSet);

%% (7) Plot the STITCH graph using precomputed XY coordinates
XY = import_gephi_xy('gephi/Wagner2018_stitch.net');
figure; stitch_plot_graph(G, XY)

%% (8) Specify node_IDs and perform graph coarse-graining
node_IDs = [];
for j = nTimePoints:-1:1
    node_IDs = [node_IDs; string(DataSet(j).celldata.cell_IDs_names)]; 
end
[G_cg, G_cg_scaff] = stitch_coarse_grain(G, index_to_graph(G, node_IDs));
graph_to_dot(adjacency(G_cg_scaff), 'directed', 0, 'filename', 'gephi/Wagner2018_cg.dot')

%% (9) Plot coarse-grained graph
XY_cg = import_gephi_xy('gephi/Wagner2018_cg.net');
figure; stitch_plot_graph_cg(G_cg, XY_cg, [], [], 'node_size', 5)
