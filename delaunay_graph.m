function G = delaunay_graph(varargin)
%DELAUNAY_GRAPH outputs a GSPBox-compatible delaunay graph from a set of
% points in R^2 or R^3.
%
%   Usage:
%       G = delaunay_graph(varargin)
%
%   Input:
%       -   It can be a N-by-d matrix with the R^d coordinates of the N 
%       points.
%       -   Otherwise, the user has to provide, in this order, N (number of 
%       points) and d (dimension of the space {2, 3}).
%
%   Output:
%       G   : A GSPBox-compatible graph structure, where G.N is the
%             number of nodes, and G.W is the graph's weighted
%             adjacency matrix.
%
%   Example:
%       G = delaunay_graph(200, 2);
%
%   Requires:
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 06 Jan 2016

%% Parse input
assert(~isempty(varargin) && nargin < 3);
if length(varargin{1}) == 1; 
    assert(nargin == 2); 
    N = varargin{1};
    d = varargin{2};
    assert(length(d) == 1);
    assert(d == 2 || d == 3);
    points = rand(N, d);
else
    assert(nargin == 1);
    points = varargin{1};
    [N, d] = size(points); 
    assert(d == 2 || d == 3);
end

%% Delaunay triangulation
dt = delaunayTriangulation(points);
edge_list = dt.edges;

%% Assemble weight matrix
spi = edge_list(:, 1);
spj = edge_list(:, 2);
dist = sqrt(sum((points(spi,:) - points(spj, :)).^2, 2));
sigma = mean(dist)^2;
spv = exp(- (dist.^2) ./ sigma);
W = sparse(spi, spj, spv, N, N);
W = W + W';

%% Create graph structure
G = struct;
G.W = W;
G = gsp_graph_default_parameters(G);
G.coords = points;
G.type = 'delaunay';

end