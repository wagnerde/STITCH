function nSig_PCs = get_nSig_PCs(X, gene_ind, nRandTrials)
% Usage: nSig_PCs = get_nSig_PCs(X, gene_ind, nRandTrials)
%
% Estimates the number of significant principal component dimensions.
% Computes a distribution of the maximum eigenvalues obtained from PCA runs 
% performed on randomized data matrices.  This distribtution is then 
% compared to the 'true' eigenvalues from PCA on the original data matrix.
% The number of true eigenvalues that are larger than the maximum random 
% eigenvalue in >95% of random trials is reported.
%
% Inputs:
% X:                Matrix size (genes x cells) of normalized counts
%
% gene_ind:         Indices of genes to consider for PCA
%
% nRandTrials:      Number of randomizations to use to generate random matrix
%                   eigenvalues

%% SETTINGS & PATHS
show_plot = true;
addpath('scTools/generic')

%% PERFORM PCA ON STANDARDIZED EXPRESSION VALUES:
z = full(zscore(X(gene_ind,:)')');

% PCA of data matrix
[~,~,eigvals] = pca(z');

% PCA of randomized matrices
for q = 1:nRandTrials
    rng(q) % set the random number seed to q (for reproducibility)
    zRnd{q} = z;
    for k = 1:size(z,1)
        zRnd{q}(k,:) = z(k,randperm(size(z,2)));
    end

    [~,~,eigvals_rnd{q}] = pca(zRnd{q}');
end

%% GENERATE EIGENVALUE HISTOGRAM AND RANDOMIZED HISTOGRAMS:
[f, x] = hist(eigvals(eigvals>0),100); 
f = f/sum(f)/(x(2)-x(1));

for q=1:nRandTrials
    [f_rnd{q}, x_rnd] = hist(eigvals_rnd{q}(eigvals_rnd{q}>0),x); 
    f_rnd{q} = f_rnd{q}/sum(f_rnd{q})/(x_rnd(2)-x_rnd(1));
end

% Evaluate eigenvalue cut-off from randomization trials (largest eigenvalue from random matrix):
for q=1:nRandTrials
    Lplus_rnd(q) = x(find(f_rnd{q}>0,1,'last')+1);
end

%% MAKE EIGENVALUE PLOT:
if show_plot

figure

% 1. Plot true eigenvalues:
h=bar(x,f); noLegend(h); set(h,'facecolor',0.5*[1 1 1],'edgecolor',0.5*[1 1 1])

% 2. Plot randomized eigenvalues:
hold on
for q=1:nRandTrials
    rnd_i = find(f_rnd{q}>0,1,'last')+3;
    h_rnd=stairs(x_rnd(1:rnd_i),f_rnd{q}(1:rnd_i),'r','linewidth',1);
    if(q>1)
        noLegend(h_rnd)
    end
end

% Annotate plot:
set(gca,'FontSize',14)
xlabel('Eigenvalue of correlation matrix')
ylabel('Eigenvalue density')
legend({'Random'})
legend('boxoff')
subtitle({['# non-random PCs = ' num2str(sum(eigvals>max(Lplus_rnd)))]},'TopLeft',[0.02 0.01]);
axis square

end

%% CALCULATE NUMBER NON-RANDOM PCs FROM >95% OF RANDOMIZATION TRIALS
for q=1:nRandTrials
    num_nonrand_PCs_list(q) = sum(eigvals>Lplus_rnd(q));
end

% Report the number of non-random PCs found in at least 95% of
% randomization trials (ie. 5th percentile value)
nSig_PCs = quantile(num_nonrand_PCs_list,0.05);


