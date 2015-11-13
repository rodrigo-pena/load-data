function [ X ] = hde( G, m )
% HDE finds an m-dimensional high-dimensional embedding of G
%
% G = (V,E,W) is a graph, G.W is its weighted adjacency matrix, and 
% G.N = |V|, i.e., the number of nodes of G.
%
%   Input parameters:
%         G         : Structure containing graph information (see GSPBox
%                     docs). The number of nodes should be accessible at 
%                     G.N, and the weighted adjacency matrix at G.W
%         m         : Number of dimensions of the subspace. 
%                     Default: min(G.N,50);
%         
%   Output parameters:
%         X         : G.N-by-m matrix whose columns span the subspace
%
%   Example:
%         G = gsp_bunny(); % (requires GSPBox: https://lts2.epfl.ch/gsp/)
%         X = hde( G, 3 );
%
%   See also: spectral_layout_2.m
%      
%   Requires: MatlabBGL (http://dgleich.github.io/matlab-bgl/)
%
%   References:
%   [1]	Y. Koren, "Drawing Graphs by Eigenvectors," Computers and
%       Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.
%
% Author: Rodrigo Pena
% Date: 13 Nov 2015
% Testing:

%% Parse input

assert(isfield(G, 'N'))
assert(isnumeric(G.N))

assert(isfield(G, 'W'))
assert(isnumeric(G.W))

if (nargin < 2) || isempty(m); m = min(G.N, 50); end

m = round(m); % Make sure m is an integer
assert((m <= G.N), 'm should be an integer in the set {1,...,G.N}');

%% Find the m-dimensional HDE of G

% Initialization
d = Inf * ones(G.N, 1);
X = zeros(G.N, m);

% Choose pivot node at random from V
p = randi(G.N);

for i = 1:m
    % Compute the i-th coordinate using BFS (matlab_bgl must be on path)
    D = bfs(G.W, p);
    
    X(:,i) = D;
    
    d = min(d, X(:,i));
    
    % Choose next pivot
    [~, p] = max(d);
end

end
