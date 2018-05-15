function var_gene_ind = get_variable_genes(DataSet, settings)
% Usage: var_gene_ind = get_variable_genes(DataSet, settings)
%
% This function wraps a series of subroutines for automated identification 
% of highly variable genes.  
% 
% INPUTS:
% DataSet         A single record/row from the DataSet structure array
% settings        Structure containing STITCH preprocessing settings
%           
% OUTPUT:
% var_gene_ind    Row indices of highly variable genes
%

%% CODE:
% if batch_flag is present, use just the first sample to find variable genes
if isfield(DataSet,'batch_flag')
    X_norm = DataSet.X_norm(:,DataSet.batch_flag==1);
    tot_counts = DataSet.tot_counts(DataSet.batch_flag==1);
else
    X_norm = DataSet.X_norm;
    tot_counts = DataSet.tot_counts;
end

% get vscores for each gene
[vscores, FF_gene, mu_gene, ~, ~, a, b] = get_vscores(X_norm, tot_counts);

% sort all genes by vscore
[vscores_sorted, vscores_rank] = sort(vscores, 'descend');
ind_not_nan = find(~isnan(vscores_sorted), 1, 'first');
var_gene_ind = vscores_rank(ind_not_nan:end); 

% apply counts filter
counts_flag = sum(X_norm > settings.minCounts, 2) >= settings.minCells;
var_gene_ind = setdiff(var_gene_ind, find(~counts_flag), 'stable');
    
% apply topVarGenes filter 
if settings.topVarGenes > 1 % topVarGenes is expressed as nGenes
    var_gene_ind = var_gene_ind(1:settings.topVarGenes);
else % topVarGenes is expressed as a fraction
    min_vscore = prctile(vscores, 100-(100*settings.topVarGenes));
    vscore_flag = vscores > min_vscore;
    sum(vscore_flag)
    var_gene_ind = setdiff(var_gene_ind, find(~vscore_flag));
end

% apply minimum gene-gene correlation filter 
var_gene_ind = get_covarying_genes(X_norm, var_gene_ind, settings.minGeneCorr);

% exclude user-specified genes & correlated genes
if ~isempty(settings.excludeGeneNames) && ~isempty(settings.allGeneNames)
    exclude_ind = get_gene_ind(settings.excludeGeneNames, settings.allGeneNames);
    exclude_ind = get_correlated_genes(exclude_ind, X_norm, settings.excludeGeneCorr, settings.excludeIter);    
    var_gene_ind = setdiff(var_gene_ind, exclude_ind, 'stable');
end

%% PLOT:
if ~settings.plot_var_genes
    return
end

% main plot
figure('position',[400  400 400 400])
p = loglog(mu_gene, FF_gene, '.', 'MarkerEdgeColor', [0.8 0.8 0.8]); noLegend(p);
hold on

% highlight subset of points
var_gene_flag = false(length(full(mu_gene)),1); var_gene_flag(full(var_gene_ind)) = true;
p2 = loglog(mu_gene(var_gene_flag), FF_gene(var_gene_flag), '.', 'MarkerEdgeColor',[0.5 0 0]); noLegend(p2);

% theory plot
hold on
x_min = 0.5*min(mu_gene(mu_gene>0));
x_max = 2*max(mu_gene);
xTh = x_min * 10.^(log10(x_max/x_min)*(0:0.01:1));
yTh = (1 + a)*(1+b) + b.*xTh;
plot(xTh,yTh,'k--','linewidth',1)

% formatting
set(gca,'fontsize',12)
xlabel('Mean Reads per Cell')
ylabel('Gene Fano Factor')
set(gca,'xlim',[min(mu_gene)*0.2 max(mu_gene)*5],'xtick',10.^(-5:5))
set(gca,'ylim',[min(FF_gene)*0.2 max(FF_gene)*5],'ytick',10.^(-5:5))
title(['Variable Gene Selection: ' DataSet.name])

