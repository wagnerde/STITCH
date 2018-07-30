function [v_scores, FF_gene, mu_gene] = get_vscores_legacy(X_use, CV_eff, CV_input)
% Usage: [v_scores, FF_gene, mu_gene] = get_vscores_legacy(X_use, CV_eff, CV_input)
% 
%  Code to estimate the v-statistic of genes from single cell RNA
%  sequencing data, following the Theory Supplement in Klein et al. (Cell,
%  2015). Coding credits to Naren Tallapragada and Allon Klein.
%
% BACKGROUND:
% The v-statistic for each gene is a corrected Fano Factor accounting for
% noise in method efficiency, and variation in cell size. Genes with highest
% v-scores are the most biologically variable. (The Fano Factor fails
% because high Fano Factor could result from variation in cell size or
% technical noise, particularly for very highly expressed genes). The
% v-stastistic is defined as:
%
%                       V = F / (A + B*mu)
%
% where mu is the mean expression of UMI-filtered data after total-count
% normalization; F=Var/mu is the Fano Factor; A = 1 + (CV_eff)^2; 
% and B = [CV_(1/N)]^2. CV_eff is the noise in the efficiency of transcript
% capture between single cells, and CV_(1/N) is the variation in total
% number of poly-A target mRNA molecules per cell in the sample. See Eq.
% (S13) in (Klein et al., Cell 2015)
%
% INPUTS:
% X_use:        Count-normalized data matrix (genes x cells).
%
% CV_eff:       Estimated CV of the efficiency. 
%
% CV_input:     Estimated CV of total mRNA input per cell, e.g. cell size
%               (Formally, CV_input = CV_(1/N), where N=total input mRNA
%
% OUTPUTS:
% v_scores:     Vector length p giving the v-statistic for each gene.
%
% FF_gene:      Vector length p giving the fano factor for each gene.
%
% mu_gene:      Vector length p giving the mean for each gene.
%

%% Calculate v-statistic
mu_gene = mean(X_use,2);
FF_gene = var(X_use,0,2)./mu_gene;

v_scores = FF_gene ./ (CV_input^2.*mu_gene + 1+CV_eff^2);
v_scores(mu_gene==0) = NaN;

