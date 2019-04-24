function EdgeL=adj2gephilab(filename,ADJ,parameters)
% Convert ana adjacency matrix of a graph to 2 spreadhseets csv files
% one for the edge table and the other for node table.
% The files _node.csv and _edge.csv have to be open 
%in Gephi via Data Laboratory.
% INPUTS:
%           filename: string for the prefix name of the two files .csv
%           ADJ: the adjacency matrix
%           parameters: vector as properties of the node to use as
%                              attributes of them.
% OUTPUTS:
%            two csv spreadsheet files: 
%                       filename_node.csv
%                       filename_edge.csv
%             EdgeL = it returns the edge list corresponing to the
%             adjacency matrix. (it can be saved to be open in Gephi too)
%
% The two files must be open in Gephi via the Data Laboratory

nodecsv=[filename,'_node.csv'];
edgecsv=[filename,'_edge.csv'];
n=size(ADJ,1); % square adjacency matrix

if nargin<3 
    parameters=ones(n,1);% all nodes have the same attributes
end
ps=parameters;
%% Node Table:
% header for node csv
fidN = fopen(nodecsv,'w','native','UTF-8'); % best format for gephi
fprintf(fidN,'%s\n','Id;Label;Attribute');
% 
for i=2:n+1
    %fprintf(fidN,'%s\n',[ num2str(i-1) ';"Node ' num2str(i-1) '"'...
    fprintf(fidN,'%s\n',[ num2str(i-1) ';"' num2str(i-1) '"'...
        ';' num2str(ps(i-1))]);
end
fclose(fidN);

%% Edge Table
EdgeL=conv_EdgeList(ADJ);

S=EdgeL(:,1); % sources
T=EdgeL(:,2); % targets
W = EdgeL(:,3); % weights

fidE = fopen(edgecsv,'w','native','UTF-8');
% header for edge csv
fprintf(fidE,'%s\n','Source;Target;Label;Weight');

for i=2:length(S)+1,
      fprintf(fidN,'%s\n',[ num2str(S(i-1)) ';' num2str(T(i-1)) ';'...
          '"Edge from ' num2str(S(i-1)) ' to ' num2str(T(i-1)) '"' ';'...
          num2str(W(i-1))]);
end
 fclose(fidE);

%% Aux function
function EdgeL=conv_EdgeList(adj)
% convert adj matrix to edge list
n=size(adj,1); % number of nodes
edges=find(adj>0); % indices of all edges
n_e=length(edges);
EdgeL=zeros(n_e,3);
for e=1:n_e
  [i,j]=ind2sub([n,n],edges(e)); % node indices of edge e  
  EdgeL(e,:)=[i j adj(i,j)];
end
