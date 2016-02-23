function [G_cell, node_cell] = connected_subgraphs(G)
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
%       G_cell      : A cell with the connected subgraphs of G
%       node_cell   : A cell with the indices of the nodes of each 
%                     connected subgraph referenced to the original graph G
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
node_cell = cell(1);

%% Find connected componnents
[n_conn_components, labels] = graphconncomp(G.W);
node_count = zeros(1, n_conn_components);

%% Assemble each subgraph structure
for i = 1:n_conn_components
    nodes = find(labels == i);
    if ~isempty(nodes)
        node_cell{i} = nodes;
        node_count(i) = length(nodes);
        G_cell{i}.W = G.W(nodes,nodes);
        if isfield(G, 'type'); G_cell{i}.type = G.type; end
        if isfield(G, 'coords'); G_cell{i}.coords = G.coords(nodes, :); end
        if isfield(G, 'lap_type'); G_cell{i}.lap_type = G.lap_type; end
        G_cell{i} = gsp_graph_default_parameters(G_cell{i});
    end
end

%% Sort the subgraphs by node count
[~, ind] = sort(node_count, 2, 'descend');
G_cell_tmp = G_cell;
node_cell_tmp = node_cell;
for i = 1:n_conn_components
    G_cell{i} = G_cell_tmp{ind(i)};
    node_cell{i} = node_cell_tmp{ind(i)};
end

end