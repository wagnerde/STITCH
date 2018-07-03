function gene_ind = get_gene_ind(query_gene_list, big_gene_list)
% Usage: gene_ind = get_gene_ind(query_gene_list, big_gene_list)
%
% Matches query gene names to a large list, returns indices for each query
%
%
%% CODE:

gene_ind = [];    
for k = 1:length(query_gene_list)
    gene_ind{k} = find(strcmp(query_gene_list{k}, big_gene_list));
end
 
if iscell(gene_ind)
    gene_ind = cell2mat(gene_ind);
end
        
