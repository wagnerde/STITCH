function var_gene_ind = get_variable_genes(DataSet, varargin)
% Usage: var_gene_ind = get_variable_genes(DataSet, varargin)
%
% This function wraps a series of subroutines for identifying and filtering
% a set of highly variable genes associated with a single-cell counts matrix.
% Includes the following steps:
% (1) Rank all genes according to a v-score statistic (a corrected Fano
%     factor).
% (2) Filter genes based on minimum counts thresholds per cell.
% (3) Filter genes based on a minimum correlation to any other gene
% (4) Exclude user-specified genes
% (5) Export a Fano-factor plot highlighting the selected genes.
%
% Inputs:
% DataSet  
%           Structure array with multiple fields and one or more records/rows 
%           for independent samples/timepoints.  
%           This function processes a single DataSet record at a time. 
%           Fields:
%           DataSet.name        A unique sample identifier string. (required) 
%           DataSet.ind         A numeric sample index (e.g. order in a
%                               time series. 
%           DataSet.X           A single-cell UMI counts matrix; (required)   
%                               rows=transcripts, columns=cells.  
%           DataSet.X_norm      A normalized UMI counts matrix. (required)   
%           DataSet.tot_counts  A vector of pre-normalized UMI counts totals
%                               for each cell. (required)           
%           DataSet.batch_flag  is a numeric index of individual sample 
%                               batches within a timepoint. 
%                               If present, gene normalizations are 
%                               performed within each batch. (optional)           
%
% Optional input name/value pairs: 
% 'topVarGenes'
%           Initial number of top variable genes to consider (default=2000)
%           Can alternatively be expressed as a fraction 
%           (e.g., 0.05 = top 5% most variable genes)
%
% 'CV_eff'
%           Optional parameter for identifying highly variable genes, passed
%           to the get_vscores_legacy function.
%           CV_eff is the noise in the efficiency of transcript capture 
%           between single cells. If left, blank, this parameter will be 
%           estimated automatically from the data (default=[])  
%
% 'CV_input'
%           Optional parameter for identifying highly variable genes, passed
%           to the get_vscores_legacy function.
%           CV_input is the variation in total number of poly-A target mRNA 
%           molecules per cell in the sample. If left, blank, this 
%           parameter will be estimated automatically from the data 
%           (default=[])  
%
% 'minCounts'
%           Filter variable genes based on minimum number of counts
%           (default=3)
%
% 'minCells'
%           Minimum number of cells that must exceed minCounts (default=1)
%
% 'minGeneCorr'
%           Filter variable genes based on minimum expression correlation
%           to other genes (default=0.2)
%
% 'excludeGeneNames'
%           Cell array of strings of gene names to exclude from the 
%           variable gene list, e.g. cell cycle markers (default={}).
%
% 'allGeneNames'
%           Cell array of strings of all gene names. Only required if 
%           'excludeGeneNames' is invoked (default={}).
%           
% 'excludeGeneCorr'
%           Expand the excluded gene list to include correlated genes based
%           on this threshold (default=0.4)
%
% 'excludeIter'
%           Number of iterations for adding correlated genes to 
%           excludeGeneList (default=2)
%
% 'show_plot'
%           Display plot of gene Fano Factor vs Mean with selected genes
%           highlighted (default=false)
%
% Output:
% var_gene_ind    Row indices of highly variable genes
%

%% PATHS
addpath('scTools/generic')

%% PARAMETER SETTINGS
% Set defaults
def.topVarGenes = 2000;
def.CV_eff = [];
def.CV_input = [];
def.minCounts = 3;
def.minCells = 1;
def.minGeneCorr = 0.2;
def.excludeGeneNames = {};
def.allGeneNames = {};
def.excludeGeneCorr = 0.4;
def.excludeIter = 2;
def.show_plot = true;

% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'get_variable_genes';
parserObj.StructExpand = false; 
parserObj.addOptional('topVarGenes',def.topVarGenes);
parserObj.addOptional('CV_eff',def.CV_eff);
parserObj.addOptional('CV_input',def.CV_input);
parserObj.addOptional('minCounts',def.minCounts);
parserObj.addOptional('minCells',def.minCells);
parserObj.addOptional('minGeneCorr',def.minGeneCorr);
parserObj.addOptional('excludeGeneNames',def.excludeGeneNames);
parserObj.addOptional('allGeneNames',def.allGeneNames);
parserObj.addOptional('excludeGeneCorr',def.excludeGeneCorr);
parserObj.addOptional('excludeIter',def.excludeIter);
parserObj.addOptional('show_plot',def.show_plot);

% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% CODE:
% if batch_flag is present, use just the first batch to find variable genes
if isfield(DataSet,'batch_flag')
    [X_norm, tot_counts] = get_normalized_counts(DataSet.X(:,DataSet.batch_flag==1));
else
    X_norm = DataSet.X_norm;
    tot_counts = DataSet.tot_counts;
end

% get vscores for each gene
% if CV_eff and CV_input are not provided, calculate automatically
if isempty(settings.CV_eff) && isempty(settings.CV_input)
    [vscores, FF_gene, mu_gene, CV_eff, CV_input] = get_vscores(X_norm, tot_counts);
else % otherwise, calculate vscores from input parameters
    CV_eff = settings.CV_eff;
    CV_input = settings.CV_input;
    [vscores, FF_gene, mu_gene] = get_vscores_legacy(X_norm, CV_eff, CV_input);    
end

% sort all genes by vscore
[vscores_sorted, vscores_rank] = sort(vscores, 'descend');
ind_not_nan = find(~isnan(vscores_sorted), 1, 'first');
var_gene_ind = vscores_rank(ind_not_nan:end); 

% apply counts filter
counts_flag = sum(X_norm >= settings.minCounts, 2) >= settings.minCells;
var_gene_ind = setdiff(var_gene_ind, find(~counts_flag), 'stable');
    
% apply topVarGenes filter 
if settings.topVarGenes > 1 % topVarGenes is expressed as nGenes
    var_gene_ind = var_gene_ind(1:settings.topVarGenes);
else % topVarGenes is expressed as a fraction
    min_vscore = prctile(vscores, 100-(100*settings.topVarGenes));
    vscore_flag = vscores > min_vscore;
    sum(full(vscore_flag));
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

%% FF PLOT:
if ~settings.show_plot
    return
end

% main plot
figure('position',[400  400 400 400])
p = loglog(mu_gene, FF_gene, '.', 'MarkerEdgeColor', [0.8 0.8 0.8]); noLegend(p);
hold on

% highlight subset of points
var_gene_flag = false(length(full(mu_gene)),1); var_gene_flag(full(var_gene_ind)) = true;
p2 = loglog(mu_gene(var_gene_flag), FF_gene(var_gene_flag), '.', 'MarkerEdgeColor',[0 0 0]); noLegend(p2);

% theory plot
hold on
x_min = 0.5*min(mu_gene(mu_gene>0));
x_max = 2*max(mu_gene);
xTh = x_min * 10.^(log10(x_max/x_min)*(0:0.01:1));
%yTh = (1 + a)*(1+b) + b.*xTh;
yTh = (1 + CV_eff^2 + CV_input^2.*xTh);
plot(xTh,yTh,'k-','linewidth',1)
legend({'CV^2= (1 + CV_e_f_f^2)/<n> + CV_i_n_p_u_t^2'})
legend('location','northwest'), legend('boxoff')
annotation('textbox',[0.15 0.84 0.3 0.01],'String',['CV_e_f_f = ', num2str(CV_eff,3)],'EdgeColor','none');
annotation('textbox',[0.15 0.80 0.3 0.01],'String',['CV_i_n_p_u_t = ', num2str(CV_input,3)],'EdgeColor','none');

% formatting
set(gca,'fontsize',12)
xlabel('Mean Reads per Cell')
ylabel('Gene Fano Factor')
set(gca,'xlim',[min(mu_gene)*0.2 max(mu_gene)*5],'xtick',10.^(-5:5))
set(gca,'ylim',[min(FF_gene)*0.2 max(FF_gene)*5],'ytick',10.^(-5:5))
title(['Gene V-scores: ' DataSet.name], 'Interpreter', 'none')

