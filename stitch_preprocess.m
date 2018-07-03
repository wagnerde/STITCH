function DataSet = stitch_preprocess(DataSet, varargin)
% Usage: DataSet = stitch_preprocess(DataSet, varargin)
%
% Preprocessing function for the STITCH pipeline.  Inputs a DataSet object 
% and performs total counts normalization and identifies highly variable
% genes.
%
% INPUTS:
% DataSet  
%           Structure array with multiple fields, and one record for each
%           timepoint. Timepoints should be ordered sequentially.
%           DataSet.ind         is a numeric index of timepoint. (required)
%           DataSet.name        is a unique string identifier. (required) 
%           DataSet.X           is a single-cell counts matrix; (required)   
%                               rows=transcripts, columns=cells.  
%           DataSet.gene_ind    is an array of gene row indices, e.g. 
%                               highly-variable genes. (optional). if
%                               provided, automated detection of highly
%                               variable genes will be skipped.
%           DataSet.batch_flag  is a numeric index of individual sample 
%                               batches within a timepoint. 
%                               if present, gene normalizations are 
%                               performed within each batch. (optional)
%                               
% Optional input name/value pairs: 
% 'topVarGenes'
%           Initial number of top variable genes to consider (default=2000)
%           Can alternatively be expressed as a fraction 
%           (e.g., 0.05 = top 5% most variable genes)
%
% 'CV_eff'
%           Optional legacy parameter for calculating gene v-scores.  
%           CV_eff is the noise in the efficiency of transcript capture 
%           between single cells. If left, blank, this parameter will be 
%           estimated automatically from the data (default=[])  
%
% 'CV_input'
%           Optional legacy parameter for calculating gene v-scores.  
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
%           on this threshold (default=0.2)
%
% 'excludeIter'
%           Number of iterations for adding correlated genes to 
%           excludeGeneList (default=2)
%
% 'plot_var_genes'
%           Display plot of gene Fano Factor vs Mean with selected genes
%           highlighted (default=false)
%

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
def.plot_var_genes = true;

% Create parser object
parserObj = inputParser;
parserObj.FunctionName = 'stitch';
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
parserObj.addOptional('plot_var_genes',def.plot_var_genes);

% Parse input options
parserObj.parse(varargin{:});
settings = parserObj.Results;

%% CODE:

% Perform total counts normalization 
nTimePoints = length(DataSet);
if (~isfield(DataSet, 'X_norm') || ~isfield(DataSet, 'tot_counts'))
    disp('Normalizing data...')
    for j = 1:nTimePoints
        [DataSet(j).X_norm, DataSet(j).tot_counts] = get_normalized_counts(DataSet(j).X);
    end
    disp('Done normalizing')
else
    disp('Normalized counts matrices already provided, skipping.')
end

% Identify variable genes
if ~isfield(DataSet, 'gene_ind')
    disp('Getting variable genes...')
    for j = 1:nTimePoints
        DataSet(j).gene_ind = get_variable_genes(DataSet(j), settings);
    end
else
    disp('Variable gene indices already provided, skipping.')
end


