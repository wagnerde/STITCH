function DataSet = stitch_preprocess(DataSet, gene_names)
% Usage: DataSet = stitch_preprocess(DataSet, gene_names)
%
% Performs data preprocessing using parameters from Wagner et al 2018.
%

%% CODE:

% Add path to scTools
addpath('scTools')

% Perform total counts normalization 
nTimePoints = length(DataSet);
disp('Normalizing data...')
for j = 1:nTimePoints
    [DataSet(j).X_norm, DataSet(j).tot_counts] = get_normalized_counts(DataSet(j).X);
end
disp('Done normalizing')

% Identify variable genes
CV_eff =   [0.60 0.40 0.40 0.40 0.40 0.40 0.40];
CV_input = [0.20 0.15 0.15 0.20 0.22 0.22 0.22];
exclude_genes = {'cdk1','mcm2', 'mcm7', 'rrm2', 'cenpa', 'cdc6', 'ccnf', 'cdca4', ...
                 'ccnd1', 'kif4','hmgb1b', 'hmgb3a', 'hspd1', 'hspa9', 'rplp0', ...
                 'hnrnpaba', 'rps2', 'rps12', 'rpl12', 'rps13', 'rps14', 'rps15a', ...
                 'rpl10', 'rps3a', 'rpl31', 'rpl37', 'rps6', 'rpl9', 'rpl11', ...
                 'rpl34', 'rpl13', 'rpl36a', 'rpl26', 'rps8a', 'rpl21', 'rps27.1',...
                 'rpl27a', 'cirbpb'};

disp('Getting variable genes...')
for j = 1:nTimePoints
    DataSet(j).gene_ind = get_variable_genes(DataSet(j), 'CV_eff', CV_eff(j) , 'CV_input', CV_input(j), 'excludeGeneNames', exclude_genes, 'allGeneNames', gene_names);
end



