function tSNE = get_simple_tSNE(X_use, gene_ind, varargin)
% Usage: tSNE = get_simple_tSNE(X_use, gene_ind, varargin)
%
% This function performs "t-distributed Stochastic Neighbor Embedding" 
% Paper:
%    Laurens van der Maaten, Geoffrey Hinton; 9(Nov):2579--2605, 2008.
%    http://www.cs.toronto.edu/~fritz/absps/tsnefinal.pdf
% Code available from:
%    https://lvdmaaten.github.io/tsne/
%
%
% INPUTS:
% X_use:        (genes x cells) matrix of counts
%
% gene_ind:    indices of genes to use for the tSNE map
% 
% Optional input name/value pairs:
% 'use_zscore':             standardize each gene individually before running 
%                       tSNE (default=true). 
%
% 'tSNE_dims':          provide optional number of dimensions to the tSNE
%                       algorithm. tSNE performs PCA on inputs and keeps
%                       the first tSNE_dims. Default = 30
% 
% 'Perplexity':         provide optional perplexity argument to tSNE.
%                       Default = 30.
%
% OUTPUTS:
% c_tSNE:       data structure with fields c_tSNE.d2, c_tSNE.d3, c_tSNE.d4
%               that are the respective tSNE maps in 2,3,4 dimensions. The 
%               n-d map is a n x L matrix, where L is the number of cells.
%
% g_tSNE:       as c_tSNE, but for genes evaluated for 2d and 3d
%

%% ADD tSNE FUNCTIONS TO PATH 
addpath('scTools/tSNE')

%% INPUTS
% Set defaults
def.use_zscore = true;
def.sampleName= '';
def.tSNE_dims = 30;
def.perplexity = 30;
def.iterations = 1000;

% Create parser object:
parserObj = inputParser;
parserObj.FunctionName = 'tSNE_cells_genes_DEW';
parserObj.StructExpand = false; % Don't treat a structure array as multiple inputs
parserObj.addOptional('sampleName',def.sampleName);
parserObj.addOptional('use_zscore',def.use_zscore);
parserObj.addOptional('tSNE_dims',def.tSNE_dims);
parserObj.addOptional('perplexity',def.perplexity);
parserObj.addOptional('iterations',def.iterations);

% Parse input options:
parserObj.parse(varargin{:});
opt = parserObj.Results;

% Currently no error checking for input options

%% GENERATE T-DISTRIBUTED SNE MAPS FOR GENES AND CELLS
% Optionally standardize gene expression (0 mean, std=1) before analysis
% (Default = standardize)
if(opt.use_zscore)
    z = full(zscore(X_use(gene_ind,:)')');
else
    z = full(X_use(gene_ind,:));
end
    
% To find embedding for genes:
%tSNE.g.d2 = tsne(z, [],2, opt.tSNE_dims, opt.perplexity, opt.iterations);   % in two dimensions
%tSNE.g.d3 = tsne(z, [],3, opt.tSNE_dims, opt.perplexity, opt.iterations);   % in three

% To get cells embedded in 2,3,4 dimensions:
tSNE.c.d2 = tsne(z', [], 2, opt.tSNE_dims, opt.perplexity, opt.iterations);
%tSNE.c.d3 = tsne(z', [], 3, opt.tSNE_dims, opt.perplexity, opt.iterations);
%tSNE.c.d4 = tsne(z', [], 4, opt.tSNE_dims, opt.perplexity, opt.iterations);


