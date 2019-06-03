function [v_scores, CV_eff, CV_input, ax_CV, ax_FF] = get_single_cell_Vscores(X_norm, tot_counts, varargin)
% [v_scores CV_eff CV_input] = get_single_cell_Vstatistic(X_norm, tot_counts, varargin)
% 
%  Code to estimate the v-statistic of genes from single cell RNA
%  sequencing data, following the Theory Supplement in Klein et al. (Cell,
%  2015).
%
% BACKGROUND:
% The v-statistic for each gene is a corrected Fano Factor accounting for
% noise in method efficiency, and variation in cell size. Gees with highest
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
% This subroutine estimates CV_(1/N) from the data. CV_(eff) is either
% independently estimated from the data, or else it is calculated from the
% relationship: 1+[CV_(eff)]^2 = (1+CV_M^2)(1+[CV_(1/N)]^2), where CV_M is
% the CV of the total number of counts detected per cell across all genes.
% CV_M is directly calculated from the data. If the user chooses to
% estimate CV_eff directly from the data, it is done as follows: the
% algorithm searches for a value of CV_eff such that the "baseline" V-score
% will be as close as possible to 1. This is done by looking for the
% maximum in the histogram of the log(Fano Factor), and setting
% (1+CV_eff^2) = F_max, where F_max is the mode of the log(Fano Factor)
% distribution.
%
% The default is to estimate CV_(1/N) and CV_eff, but an option allows
% fitting CV_(1/N) only and inferring CV_eff from CV_M.
%
% Because CV_(1/N) and CV_eff are estimates, the v-score is not uniquely
% defined and it's precise value will change depending on how the
% estimation is carried out. Several options allowing tweaking the fitting
% algorithm.
% 
% INPUTS:
% X_norm            : Matrix of total-count normalized single cell gene
%                     expression data with size p x q, where p=(# genes),
%                     q=(# cells).
% 
% tot_counts        : Vector length q with the total number of UMI-filtered
%                     mapped transcripts per cell across all genes, before
%                     normalization.
%
% varargin          : See below for optional name/value pairs to tweak
%                     estimation of CV_(1/N) and resulting v-scores, and
%                     for optional visualization.
%
% OUTPUTS:
% v_scores          : vector length p giving the v-statistic for each gene.
%
% CV_eff            : Estimated value of CV_eff from the fit to the data.
%
% CV_input          : Estimated value of CV_(1/N) from the fit to the data.
%
% figCV             : Handle to figure, if requested
%
% figFF             : Handle to figure, if requested
%
% OPTIONAL NAME-VALUE PAIRS:
% 'fit_CVeff'       : Logical. If true, CVeff is treated as an independent
%                     fitting parameter. If false, only CV_(1/N) is fitted
%                     with  CV_eff^2 = (1+(CV_M)^2)(1+(CV_(1/N))^2)-1.
%                     Default = false.
%
% 'min_mean'        : Scalar z, double. Only include genes with average
%                     expression greater than z. Default z=0.
%
% 'fit_percentile'  : Scalar f, double. Fit baseline curve (V=1) to this
%                     specified running percentile value as function of
%                     mean. Default is f=0.1
%
% 'error_wt'        : Scalar wt. The exponent of fitting function. A value
%                     of wt=2 is least-squares (L2). wt=1 is robust (L1).
%                     Default wt=1.
%
% 'show_plot'       : Optionally plot CV vs mean with genes color-coded by
%                     v-score.
%
% 'saturation_v'    : If plotting, this sets the value of the v-score at
%                     which the color bar saturates. Default = 10. Set to
%                     very large or to Inf to remove saturation.

%% PATHS
addpath('scTools/generic')

%% DEFAULTS 
% Set defaults
def.fit_CVeff       = true;
def.min_mean        = 0; % Exclude genes with average expression of this value or lower
def.fit_percentile  = 0.33; % Fit to this percentile of CV-vs-mean values
def.error_wt        = 1; % Exponent of fitting function. Value of 2 is least-squares (L2). Value of 1 is robust (L1).
def.show_plot       = false; % Plot CV vs mean plot with baseline curve and color-coded by v score
def.saturation_v    = 10; % Visualization: value of V above which colorbar saturates

% Create parser object:
parserObj = inputParser;
parserObj.FunctionName = 'get_single_cell_Vscores';
parserObj.addOptional('fit_CVeff',def.fit_CVeff);
parserObj.addOptional('min_mean',def.min_mean);
parserObj.addOptional('error_wt',def.error_wt);
parserObj.addOptional('fit_percentile',def.fit_percentile);
parserObj.addOptional('show_plot',def.show_plot);
parserObj.addOptional('saturation_v',def.saturation_v);
    
% Parse input options:
parserObj.parse(varargin{:});
opt = parserObj.Results;


%% PREPARE DATA FOR FITTING
mu_gene   = mean(X_norm,2);
FF_gene   = ( std(X_norm,0,2).^2 ) ./ mu_gene;

% Perform fit on log-transformed data:
data_x    = log(mu_gene);
data_y    = log(FF_gene./mean(X_norm,2));

% Exclude genes with 0 counts (or fewer than user-defined "min_mean")
ind       = find(mean(X_norm,2)<=opt.min_mean);
data_x_fit= data_x; data_x_fit(ind) = [];
data_y_fit= data_y; data_y_fit(ind) = [];

% GENERATE A RUNNING PERCENTILE OF THE FITTING DATA:
nBins = 50;
[x, y] = runningquantile(data_x_fit,data_y_fit,...
    opt.fit_percentile, nBins);

% ELIMINATE ANY NANS PRODUCED FROM HAVING BINS THAT ARE EMPTY:
x(isnan(y)) = [];
y(isnan(y)) = [];

% Data is now ready for fitting.

%% DEFINE FITTING FUNCTION to evaluate with values from parameter fitting
% For "baseline", genes are non-variable and have v-statistic V=1. For
% these genes, CV^2 = A/mu + b, where A=(1+a)(1+b), and a=CV_M^2,
% b=CV_(1/N)^2.

% For fitting CV_(1/N) only:
%f = @(x,a,b)((1+a)*(1+b)./x + b); 
fLog = @(x,a,b)(log((1+b)*(1+a)*exp(-x) + b));

% For fitting CV_(1/N) and CV_eff:
%g = @(x,c,b)(c./x + b);
gLog = @(x,c,b)(log(c*exp(-x) + b));

%% CALCULATE a=(CV_M)^2 from unnormalized UMI counts
if(~opt.fit_CVeff)
    a = (std(tot_counts)/mean(tot_counts))^2;
else
    [FFhist FFhist_vals] = histcounts(log10(FF_gene(mu_gene>0)),200);
    c = exp(FFhist_vals(FFhist==max(FFhist))); % z = 1/( (1+a)(1+b) )
    c = c(1); % In case of ambiguous maximum
    c = max([1 c]); % Do not allow c to fall below 1.
end
%% PERFORM FITTING
wt = opt.error_wt;
if(~opt.fit_CVeff)
    % FIT b=[CV_(1/N)]^2
    errFun = @(b)(sum(abs(fLog(x,a,b)-y).^wt));
    b0 = 0.1;           % Make a starting guess at the solution
    options=optimset('Display','off','MaxFunEvals',1000,'Algorithm','interior-point');   % Option to display output
    b = fminsearch(errFun,b0,options);
else
    % FIT c=[CV_eff]^2
    errFun = @(b)(sum(abs(gLog(x,c,b)-y).^wt));
    b0 = [0.1];           % Make a starting guess at the solution
    options=optimset('Display','off','MaxFunEvals',1000,'Algorithm','interior-point');   % Option to display output
    b = fminsearch(errFun,b0,options);
    a = c/(1+b) - 1;
end
%b=0.6^2;    
%% Compute v-score test statistic for all genes
v_scores    = FF_gene ./ ((1+a)*(1+b) + b * mu_gene);

%% OUTPUT
CV_eff = sqrt((1+a)*(1+b) - 1);
CV_input = sqrt(b);

%% OPTIONAL PLOT
if (~opt.show_plot)
    return
end


figure('position',[360   286   465   412],'Name','Gene CV vs mean and V statistic')
% y values
CV_gene = sqrt(FF_gene./mu_gene);
% Color values:
plt_v = v_scores;
plt_v(plt_v>opt.saturation_v) = opt.saturation_v;

% Plot:
p=scatter(mu_gene, CV_gene,2, plt_v); noLegend(p);
colormap(jet(100))
set(gca,'yscale','log','xscale','log');
axis square
cbar_hndl = colorbar;
ylabel(cbar_hndl,'V statistic')
if(max(v_scores)>opt.saturation_v)
    cbyl = get(cbar_hndl,'yticklabel');
    cbyl{end} = ['>' num2str(opt.saturation_v)];
    set(cbar_hndl,'yticklabel',cbyl)
end

% Theory plots:
hold on
x_min = 0.5*min(mu_gene(mu_gene>0));
x_max = 2*max(mu_gene);
xTh = x_min * 10.^(log10(x_max/x_min)*[0:0.01:1]);
yTh = sqrt((1 + a)*(1+b)./xTh + b);
plot(xTh,1./xTh.^0.5,'-','color',[1 0.5 0.5],'linewidth',1)
plot(xTh,yTh,'r-')
set(gca,'fontsize',16)
xlabel('avg. reads per cell')
ylabel('Coefficient of Variation')
set(gca,'xlim',[1e-2 x_max],'xtick',10.^[-4:4])
set(gca,'ylim',[0.9*min(CV_gene(mu_gene>0)) 1.1*max(CV_gene(mu_gene>0))],'ytick',10.^[-1:3])

legend({'CV^2=1/<n>',['CV^2= (1+' num2str(CV_eff,3) '^2)/<n> + ' num2str(CV_input,2) '^2']})
legend('location','northeast'), legend('boxoff')

ax_CV = gca;

%% Fano factor plot
figure('position',[360   286   465   412],'Name','Gene FF vs mean and V statistic')
% Plot:
p=scatter(mu_gene, FF_gene,2, plt_v); noLegend(p);
colormap(jet(100))
set(gca,'yscale','log','xscale','log');
axis square
cbar_hndl = colorbar;
ylabel(cbar_hndl,'V statistic')
if(max(v_scores)>opt.saturation_v)
    cbyl = get(cbar_hndl,'yticklabel');
    cbyl{end} = ['>' num2str(opt.saturation_v)];
    set(cbar_hndl,'yticklabel',cbyl)
end

% Theory plots:
hold on
x_min = 0.5*min(mu_gene(mu_gene>0));
x_max = 2*max(mu_gene);
xTh = x_min * 10.^(log10(x_max/x_min)*[0:0.01:1]);
yTh = (1 + a)*(1+b) + b.*xTh;
plot(xTh,yTh,'r-')
set(gca,'fontsize',16)
xlabel('avg. reads per cell')
ylabel('Fano Factor')
set(gca,'xlim',[1e-2 x_max],'xtick',10.^[-4:4])
set(gca,'ylim',[0.9*min(FF_gene(mu_gene>0)) 1.1*max(FF_gene(mu_gene>0))],'ytick',10.^[-1:3])

legend({['CV^2= (1+' num2str(CV_eff,3) '^2)/<n> + ' num2str(CV_input,2) '^2']})
legend('location','northeast'), legend('boxoff')

ax_FF = gca;

end

%% SUPPORT FUNCTION
function [xOut, yOut] = runningquantile(x,y,p,nBins)
% [xOut, yOut] = runningquantile(x,y,p,nBins)
% Calculates a running p-th quantile of y over x.
% (x, y) = 2 same-length vectors of input data
% p = between 0-1, the quantile
% nBins = number bins to divide x range into.
%
% Output:
% xOut: x values at regular intervals
% yOut: average y value for each value of xOut
% points in bin)

[x ind]=sort(x);
y = y(ind);

dx = (x(end)-x(1))/(nBins);
xOut = (x(1)+dx/2):dx:(x(end)-dx/2);

yOut = zeros(size(xOut));
for i=1:length(xOut)
    ind = find(x(:)>=xOut(i)-dx/2 & x(:)<xOut(i)+dx/2);
    if(~isempty(ind))
        yOut(i) = quantile(y(ind),p);
    else
        if(i>1)
            yOut(i) = yOut(i-1);
        else
            yOut(i) = nan;
        end
    end

end

end
