function [ G ] = higgs_social_graph(HIGGS_DIR, directed)
% Outputs a GSP Toolbox-compatible graph from the Digg 2009 dataset.
%   Inputs:
%       HIGGS_DIR:  Directory where the higgs-social-network.edgelist and
%                   the higgs-activity_time.txt files are located.
%       directed:   0 or 1(default): assemble and undirected
%                   graph, or a directed graph, respectively.
%
%   Outputs:
%       G:          A GSP Toolbox-compatible graph assembled from
%                   the Higgs-Twitter social network information.
%
%   Example:
%       home = getenv('HOME');
%       HIGGS_DIR = strcat(home, '/data/higgs-twitter');
%       [ G ] = higgs_social_graph(HIGGS_DIR);
%
% Author: Rodrigo Pena (rodrigo.pena@epfl.ch)
% Date: 30 Oct 2015

%% Process input
if nargin < 2 || isempty(directed) || (directed ~= 1 && directed ~= 0)
    directed = 1;
end

%% Read data file
fid = fopen(strcat(HIGGS_DIR, '/higgs-social_network.edgelist'));
C = textscan(fid, '%f%f');
fclose(fid);

%% Create sparse adjacency matrix for the network
spi = C{1};
spj = C{2};
N = length(spi);
spv = ones(size(spi));
A = sparse(spi,spj,spv,N,N);
A = spdiags(zeros(N,1), 0, A); % Remove self-friendships
if ~directed
    A = max(A, A');
end

%% Create graph structure
G.N = N;
G.W = A;
G.coords = [randperm(N); randperm(N)]';
G.A = logical(A);
G.type = 'unknown';
G.directed = directed;
G.d = sum(G.W,2);
G.Ne = nnz(G.W)*(directed + 0.5*~directed);
G.L = diag(G.d) - G.W;
G.lap_type = 'combinatorial';

try
    G = gsp_graph_default_plotting_parameters(G);
catch
    warning('GSPBox not found. Could not set default plotting parameters.')
end

end

