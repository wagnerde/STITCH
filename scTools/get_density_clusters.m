function clustIDs = get_density_clusters(xy)
% Usage: clustIDs = get_density_clusters(xy)
%
% This function performs "density-based clustering" on 2d point clouds.
% Paper:
%    Alex Rodriguez & Alessandro Laio: Clustering by fast search and find of density peaks,
%    Science 344, 1492 (2014); DOI: 10.1126/science.1242072.
% Code available from:
%    https://www.mathworks.com/matlabcentral/fileexchange/53922-densityclust
%       
% 

%% CODE:

% default settings
percNeigh = 0.03; 
distance_cutoff = 10;
minPoints = 10;
kernel = 'Gauss';

% add functions to path
addpath('scTools/densityClust')

% perfom density clustering
distance_matrix = squareform(pdist(xy, 'euclidean'));
[~, rho] = paraSet(distance_matrix, percNeigh, kernel);
[~, clustIDs, ~, ~] = densityClust(distance_matrix, distance_cutoff, rho, 0, minPoints);



