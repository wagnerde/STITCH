function z = get_zscore_batch(X, batch_flag)
% Usage: z = get_zscore_batch(X, batch_flag)
%
% Convert gene counts to z-scores independently for separate batches of
% cells
%
% INPUTS:    
% X            Counts matrix (rows=genes, columns=cells) for constructing
%              gene zscores
% batch_flag   A numeric index batch IDs for each cell (column) in X
%
% OUTPUT
% z            Counts matrix with rows normalized by z-score
% 

%% CODE:

batches = unique(batch_flag);   
z = zeros(size(X));
for k = 1:length(batches)
    batch_ind = batch_flag == k;
    Xb = X(:,batch_ind);
    Xb_z = zscore(Xb,0,2);
    z(:,batch_ind) = Xb_z;
end


