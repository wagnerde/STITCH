function all_node_positions = get_kNN_trajectory(G, root_nodes, term_nodes, nRandTrials, fSubSample, fKeep, exclude_nodes)
% Usage: all_node_positions = get_kNN_trajectory(G, root_nodes, term_nodes, nRandTrials, fSubSample, fKeep, exclude_nodes)
%
% This function calculates a guided linear trajectory through a kNN graph 
% in a manner inspired by Wanderlust (Bendall et al Cell 2014).  The user 
% specifies two sets of graph nodes to act as the approximate root and 
% terminus of the trajectory.  The algorithm then computes a set of
% shortest paths between each root node and each terminus node after 
% randomly subsampling a fraction of total graph edges.  This process is 
% repeated several times.  Finally, all nodes discovered by this process
% are then ordered based on mean path positions over all trials. 
%
% INPUTS:
% 'G'                 
%         A standard Matlab graph object.
%
% 'root_nodes'        
%         Indices of graph nodes to act as the root of the trajectory.
%
% 'term_nodes'        
%         Indices of nodes to act as the trajectory terminus.
%
% 'nRandTrials'       
%         Number of graph randomizations to perform.
%
% 'fSubSample'        
%         Fraction total edges to randomly remove during each trial.
%
% 'fKeep'           
%         Fraction of random trials in which a node must be 
%         discovered to retain it in the final trajectory: 
%         fKeep = 1 keeps all nodes
%         fKeep = 0.9 keeps nodes discovered in >=10 % of trials
%         fKeep = 0.5 keeps nodes discovered in >=50 % of trials
%         fKeep = 0.1 keeps nodes discovered in >=90 % of trials
%
% 'exclude_nodes'     
%         Indices of graph nodes to exclude (optional).
%
%
% OUTPUTS: 
% 'all_node_positions'
%       A numeric array of trajectory positions for each node in G.
%       Undiscovered and excluded nodes receive a value of 'NaN'
%

%% CODE:

% if specified, remove edges to any excluded nodes
if exist('exclude_nodes', 'var')
    edge_flag = ismember(str2num(cell2mat(G.Edges.EndNodes(:,1))), exclude_nodes) | ismember(str2num(cell2mat(G.Edges.EndNodes(:,2))), exclude_nodes);
    G = rmedge(G, find(edge_flag));
end
nEdges = length(G.Edges.EndNodes);
nNodes = length(G.Nodes.Name);

% find shortest paths between root and all terminal nodes
path_nodes = [];  % running list of all node discoveries
path_node_positions = [];  % running list of trajectory positions for all node discoveries
path_node_each_trial = []; % running tally of nodes discovered in each trial

for k = 1:nRandTrials % random trials
    for l = 1:length(root_nodes) % number of term1 nodes

        % remove a random subset of all graph edges
        rng(k); % sets the random number seed to k, for reproducibility
        rand_edge_ind = randperm(nEdges);
        rand_edge_ind = rand_edge_ind(1:ceil(fSubSample*nEdges));
        G_rand = rmedge(G,rand_edge_ind);

        % identify the shortest path tree connecting the root to all terminal nodes
        [TR_paths] = shortestpathtree(G_rand, term_nodes, root_nodes(l), 'OutputForm','vector'); % vector output is much faster than cell output
        path_nodes_next = find(~isnan(TR_paths))';
        sources = TR_paths(~isnan(TR_paths))';
        nPathNodes = length(path_nodes_next);

        % extract path positions for each node
        path_node_positions_next = zeros(1,nPathNodes);
        path_node_positions_next((sources == 0)) = 1;
        for j = 2:nPathNodes
            path_node_positions_next(ismember(sources,path_nodes_next(path_node_positions_next == j-1))) = j;
        end
        path_node_positions_next = path_node_positions_next / max(path_node_positions_next);

        % append path node and position vectors
        path_nodes = [path_nodes; path_nodes_next];        
        path_node_positions = [path_node_positions; path_node_positions_next'];     
    end
    
    % keep track of which nodes were discovered in this trial
    path_node_each_trial = [path_node_each_trial; unique(path_nodes)];

end

% find all unique nodes and the mean path position for each
[unique_path_nodes, ~, ic] = unique(path_nodes);
nUniquePathNodes = length(unique_path_nodes);
for j = 1:nUniquePathNodes
    unique_path_node_positions(j) = mean(path_node_positions(ic == j));
end

% keep any node that was discovered in the top fKeep fraction of random trials
keep_ind = count_unique(path_node_each_trial) > (1-fKeep)*nRandTrials;

% output trajectory positions for all nodes in the original graph 
traj_node_ind = unique_path_nodes(keep_ind);
traj_node_pos = unique_path_node_positions(keep_ind);
all_node_positions = nan(nNodes,1);
all_node_positions(traj_node_ind) = traj_node_pos;

