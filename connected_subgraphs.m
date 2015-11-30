function G_cell = connected_subgraphs(G)
%CONNECTED_SUBGRAPH outputs a cell whose entries are the connected
%subgraphs in graph G, sorted by decreasing node count.
%
%   Usage:
%       G_cell = connected_subgraphs(G)
%
%   Input:
%       G   : A Matlab structure with graph information. It must have as
%             fields:
%           G.N     : The number of nodes
%           G.W     : The graph's G.N-by-G.N weighted adjacency matrix.
%
%   Output:
%       G_cell  : A cell with the connected subgraphs of G
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   Example:
%       G_cell = connected_subgraphs(G);
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 30 Nov 2015

%% Parse input
assert(isfield(G, 'N') && isfield(G, 'W'), ...
    'G must have N (nb of nodes) and W (adj. matrix) as fields');

%% Initialization
G_cell = cell(1);

%% Dulmage-Mendelsohn permutation
[p, ~, r, ~] = dmperm(G.W); 
n_conn_components = length(r) - 1;
node_count = zeros(1, n_conn_components);

%% Assemble each subgraph structure
for i = 1:n_conn_components
    nodes = p(r(i):r(i+1)-1);
    node_count(i) = length(nodes);
    G_cell{i}.W = G.W(nodes,nodes);
    G_cell{i}.type = 'from weight';
    G_cell{i} = gsp_graph_default_parameters(G_cell{i});
end

%% Sort the subgraphs by node count
[~, ind] = sort(node_count, 2, 'descend');
G_cell_copy = G_cell;
for i = 1:n_conn_components
    G_cell{i} = G_cell_copy{ind(i)};
end

end