function L = compute_graph_laplacian(W, laptype)
% COMPUTE_GRAPH_LAPLACIAN computes the graph Laplacian from a graph's 
% weight matrix 
% 
%   Input parameters:
%         W         : Weight (or Adjacency) matrix of a graph
%         laptype   : Type of graph Laplacian to compute, 'combinatorial',
%                     or 'normalized'.
%         
%   Output parameters:
%         L         : The graph Laplacian
%
%   Example:::
%         G = gsp_swiss_roll(); (requires GSPBox)
%         L = compute_graph_laplacian(G.W);
%
%   See also: compute_graph_gradient.m
%
%   References:
%   [1]	F. R. K. Chung, ?Lectures on Spectral Graph Theory,? 2001, 
%       pp. 1-25.
%
% Author: Rodrigo Pena
% Date: 6 Nov 2015
% Testing:

%% Parse parameters
assert(size(W,1)==size(W,2), 'W is not a square matrix');

if (nargin < 2) || isempty(laptype)
    laptype = 'combinatorial';
end

assert(strcmp(laptype, 'combinatorial') || strcmp(laptype, 'normalized'), ...
    'laptype must be `combinatorial` or `normalized`');

%% Compute Laplacian
n = size(W, 1);
D = sum(W, 2);
switch laptype
    case 'combinatorial'
        L = diag(D) - W;
    case 'normalized'
        invDhalf = spdiags(1./sqrt(D + eps), 0, n, n);
        L = speye(n) - (invDhalf * W * invDhalf);
end

if ~issparse(L)
    L = sparse(L);
end

if ~issymmetric(L)
    L = (L + L') / 2;
end

end