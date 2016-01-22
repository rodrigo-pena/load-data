function LG = line_graph(G)
%LINE_GRAPH creates line graph from graph G
%
%   Usage:
%       LG = line_graph(G)
%
%   Input:
%       G   : A Matlab structure, where G.W is a N-by-N weighted
%             adjacency matrix, where N is the number of nodes in
%             graph G.
%
%   Output:
%       LG  : A Matlab structure, where LG.W is a Ne-by-Ne weighted
%             adjacency matrix, where Ne is the number of edges in
%             graph G.
%
%   Example:
%       LG = line_graph(G);
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   Reference: https://en.wikipedia.org/wiki/Line_graph
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 15 Jan 2016

%% Parse input
assert(isfield(G, 'W'), ...
    'G must have at least an adjacency matrix accessible on field G.W');
if (~isfield(G, 'N') || ~isfield(G, 'Ne')) 
    G = gsp_graph_default_parameters(G);
end

%% Compute edge list
[source, target, weight] = find(tril(G.W));

%% Assemble G's incidence matrix
spi = [1:G.Ne, 1:G.Ne]';
spj = [source; target];
spv = [-sqrt(weight); sqrt(weight)]; % So that I*I' = G.L
I = sparse(spi, spj, spv, G.Ne, G.N)';

%% Assemble the weighted adjacency matrix of the line graph
W = I'*I - sparse(diag(weight));

%% Create graph structure for LG
LG = struct('W', W);
LG = gsp_graph_default_parameters(LG);
LG.coords = (G.coords(source, :) + G.coords(target, :))/2;
LG.type = 'line_graph';

end