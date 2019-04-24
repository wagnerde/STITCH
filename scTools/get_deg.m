function DEG_RESULTS = get_deg(X, groupIDs, names_genes, varargin)
% Usage: DEG_RESULTS = get_deg(X, groupIDs, names_genes, varargin)
% 
% This function identifies differentially expressed genes for a series of 
% cell groups. In default mode, a single test is performed for each given
% groupID, comparing cells of the groupID to all other cells in the dataset.
% If 'test_groups' is specified, explicit pairs of groups are compared.
% Normalized UMI counts are converted to TPM; significance is determined by 
% a two-sided Wilcoxon rank-sum test (also known as the Mann-Whitney U test). 
% P values are adjusted by mafdr. The function reports tables of both 
% up-regulated and down-regulated genes for each test.  
% 
% INPUTS:
% 'X' 
%       A total-count-normalized transcripts x cells counts matrix 
%       (rows=transcripts; columns=cells).
%
% 'groupIDs'      
%       Numerical array of group assignments for each cell / each row in xy
%               
% 'names_genes'
%       Cell array of strings. Names for all genes (corresponding to rows in X).
%
% OPTIONAL NAME/VALUE PAIRS:
% 'max_nCells_per_test' 
%       Maximum number of cells to use for each test group. Groups that 
%       exceed this number of cells are randomly subsampled. Subsampling
%       can improve speed/performance when large numbers of cells are
%       analyzed (default=Inf).
%
% 'tpm_thresh'
%       Minimum average expression threshold for filtering/removing 
%       low-expressed transcripts (default=0).   
%
% 'fdr_thresh'
%       Significance threshold (alpha) applied to fdr-corrected p-values
%       (default=0.05).
% 
% 'log2fc_thresh'
%       Log2 fold change threshold (default=1 e.g. 2-fold change). 
% 
% 'test_groups'
%       Cell array of two sets of groupIDs, for specifying explicitly 
%       defined test groups. Alternative to default mode. When invoked, 
%       only a single test is performed.
%       e.g. {[1 2],[3]} compares cells of both groups 1 and 2 to cells of 
%       group 3 (default={}). 
%
% OUTPUTS:
% DEG_RESULTS
%       Structure array in which each row corresponds to a single DEG test.  
%       Test parameters and DEG results are reported.  Results tables
%       indclude gene names, fold changes, fdr-corrected p-values, and mean 
%       TPM expression levels for all identified DE genes. 
%

%% PARAMETER SETTINGS
% Set defaults
def.max_nCells_per_test = Inf;
def.tpm_thresh = 0;
def.fdr_thresh = 0.05;
def.log2fc_thresh = 1;
def.test_groups = {};

% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'get_DEG';
parserObj.StructExpand = false; 
parserObj.addOptional('max_nCells_per_test',def.max_nCells_per_test);
parserObj.addOptional('tpm_thresh',def.tpm_thresh);
parserObj.addOptional('fdr_thresh',def.fdr_thresh);
parserObj.addOptional('log2fc_thresh',def.log2fc_thresh);
parserObj.addOptional('test_groups',def.test_groups);

% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% CODE:

% convert UMI counts to TPM 
X_tpm = bsxfun(@rdivide,X,sum(X,1))*1e6;

% ignore any '-1' or 'NaN' assignments 
groupIDs(groupIDs==-1) = NaN; 
unique_IDs = unique(groupIDs);
unique_IDs(isnan(unique_IDs)) = []; 

% ignore genes below the specified tpm threshold
mean_tpm = mean(X_tpm,2);
exp_genes = mean_tpm > settings.tpm_thresh;
X_tpm = X_tpm(exp_genes,:);
names_genes = names_genes(exp_genes);
mean_tpm = mean_tpm(exp_genes);

% decide on the number of tests
if isempty(settings.test_groups)
    % default mode
    nTests = length(unique_IDs); 
else
    % one explicitly defined test
    nTests = 1;
end

% compare cells between test groups; compute fold changes and p-values
nGenes = size(X_tpm,1);
fch = zeros(nGenes, nTests);
pval = ones(nGenes, nTests);
fdr = ones(nGenes, nTests);

% perform tests
for j = 1:nTests
    
    % set up the two test groups
    if isempty(settings.test_groups) 
        % default mode
        disp(['Getting DEG: All vs. Group ' num2str(j) ' ...'])
        test_group1 = find(groupIDs == j);
        test_group2 = find(groupIDs ~= j); 
        DEG_RESULTS(j).test_group1 = num2str(j);
        DEG_RESULTS(j).test_group2 = 'All';   
    else
        % test groups explicitly defined
        disp(['Getting DEG: Groups ' num2str(settings.test_groups{1}) ' vs. Groups ' num2str(settings.test_groups{2}) ' ...'])
        test_group1 = find(ismember(groupIDs, settings.test_groups{1}));
        test_group2 = find(ismember(groupIDs, settings.test_groups{2}));        
        DEG_RESULTS(j).test_group1 = num2str(settings.test_groups{1});
        DEG_RESULTS(j).test_group2 = num2str(settings.test_groups{2});    
    end
    
    % if group sizes exceed max_nCells_per_test, then subsample 
    rng(802); % sets the random number seed for reproducibility
    if length(test_group1) > settings.max_nCells_per_test
        rand_ind_group1 = randperm(length(test_group1));
        test_group1 = test_group1(rand_ind_group1(1:settings.max_nCells_per_test));
    end    
    if length(test_group2) > settings.max_nCells_per_test
        rand_ind_group2 = randperm(length(test_group2));
        test_group2 = test_group2(rand_ind_group2(1:settings.max_nCells_per_test));
    end
    
    % report group sizes in the final results table
    DEG_RESULTS(j).size_group1 = num2str(length(test_group1));
    DEG_RESULTS(j).size_group2 = num2str(length(test_group2));
    
    % perform rank-sum test
    for k = 1:nGenes  
        fch(k,j) = log2((mean(X_tpm(k,test_group1))+1) / (mean(X_tpm(k,test_group2))+1));
        pval(k,j) = ranksum(X_tpm(k,test_group1), X_tpm(k,test_group2));
    end    
    fdr(:,j) = mafdr(pval(:,j));

end

disp('Done!')


% tabulate & sort results

for j = 1:nTests       
    
    % report test parameters in the final results table
    DEG_RESULTS(j).tpm_thresh = settings.tpm_thresh;
    DEG_RESULTS(j).fdr_thresh = settings.fdr_thresh;
    DEG_RESULTS(j).log2fc_thresh = settings.log2fc_thresh;
    
    % UPREGULATED GENES
    % get genes that passed thresholds
    up_this_group = find((fdr(:,j) < settings.fdr_thresh) & (fch(:,j) > settings.log2fc_thresh));    
    % sort by fold change
    [~,up_sort_ind] = sort(fch(up_this_group,j), 'descend');    
    % make a table
    up_genes_names = names_genes(up_this_group(up_sort_ind))';  
    up_genes_mean_tpm = mean_tpm(up_this_group(up_sort_ind));   
    up_genes_fch = fch(up_this_group(up_sort_ind),j);  
    up_genes_fdr = fdr(up_this_group(up_sort_ind),j);  
    DEG_RESULTS(j).UP_table = table(up_genes_names, up_genes_fch, up_genes_fdr, up_genes_mean_tpm, 'VariableNames', {'GeneName', 'Log2FC','FDR', 'avgTPM'});

    % DOWNREGULATED GENES
    % get genes that passed thresholds
    down_this_group = find((fdr(:,j) < settings.fdr_thresh) & (fch(:,j) < settings.log2fc_thresh));    
    % sort by fold change
    [~,down_sort_ind] = sort(fch(down_this_group,j), 'ascend');    
    % make a table
    down_genes_names = names_genes(down_this_group(down_sort_ind))';   
    down_genes_mean_tpm = mean_tpm(down_this_group(down_sort_ind));   
    down_genes_fch = fch(down_this_group(down_sort_ind),j);  
    down_genes_fdr = fdr(down_this_group(down_sort_ind),j);  
    DEG_RESULTS(j).DOWN_table = table(down_genes_names, down_genes_fch, down_genes_fdr, down_genes_mean_tpm, 'VariableNames', {'GeneName', 'Log2FC','FDR', 'avgTPM'});

end