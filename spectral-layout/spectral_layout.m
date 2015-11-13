function [ G, U, e ] = spectral_layout( G, indices, method, verbose )
% SPECTRAL_LAYOUT assigns new coordinates to the vertices of graph G based
% on the theory of spectral layout.
%
%   Input parameters:
%         G         : Structure containing graph information (see GSPBox
%                     docs). The graph Laplacian eigenvectors should be
%                     accessible at G.U, and the respective eigenvalues at
%                     G.e.
%         indices   : Vector of integers whose elements lie in {1,...,G.N}.
%                     If U is the matrix of eigenvectors of the graph
%                     Laplacian, the new coordinates of the vertices will
%                     be G.coords = U(:,indices).
%         method    : - 'standard': use standard eigenvectors of G.L, the
%                       graph Laplacian.
%                     - 'deg-normalized': use generalized eigenvectors of
%                       (G.L, D), where D is the diagonal degree matrix.
%         verbose   : Level of information display. Default: 0.
%
%   Output parameters:
%         G         : Structure containing graph information (see GSPBox
%                     docs), with the nodes' coordinates updated.
%         U         : The computed eigenvectors.
%         e         : The computed eigenvalues.
%
%   Example:::
%         G = gsp_bunny(); % Requires GSPBox
%         figure();
%         subplot(121)
%         gsp_plot_graph(G);
%         title('Original Coordinates')
%         indices = [2,3,4];
%         subplot(122)
%         spectral_layout(G, indices);
%
%   See also: spectral_layout_approx.m
%
%   Requires:   - GSPBox (https://lts2.epfl.ch/gsp/)
%
%   References:
%   [1]	Y. Koren, "Drawing Graphs by Eigenvectors," Computers and
%       Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.
%
% Author: Rodrigo Pena
% Date: 6 Nov 2015
% Testing: test/test_spectral_layout.m

%% Parse input

% indices
assert(length(indices)<=G.N, 'indices should not exceed G.N in length');

assert(sum(indices>G.N)==0 && sum(indices<1)==0,...
    'indices should not have elements out of the set {1,...,G.N}');

indices = round(indices); % Make sure the elements in indices are integers

% method
if (nargin < 3) || (isempty(method))
    method = 'standard';
end

switch method
    case 'standard'
        D = eye(G.N);
    case 'deg-normalized'
        D = diag(sum(G.W, 2));
    otherwise
        error('Options for method are ''standard'', or ''deg-normalized''');
end

% verbose
if (nargin < 4) || (isempty(verbose))
    verbose = 0;
end

%% Compute eigenvectors

if ~issparse(G.L) || (G.N <= 3000)
    try
        [U, e] = eig(full(G.L),full(D));
    catch
        [U, e, ~] = svd(full(D\G.L));
    end
else
    try
        [U, e] = eigs(G.L, D, max(indices), 'SM');
    catch
        [U, e, ~] = svds(D\G.L, max(indices), 0);
    end
end

e = diag(e);
[e, idx] = sort(e);
U = U(:,idx);

% Choose non-negative representation
signs = sign(U(1,:));
signs(signs==0)=1;
U = U * diag(signs);

%% Assign new coordinates
G.coords = U(:, indices);
G = gsp_graph_default_plotting_parameters(G);
minVec = min(G.coords, [], 1);
maxVec = max(G.coords, [], 1);
G.plotting.limits = reshape([minVec; maxVec], 1, []);

%% Output result
if (nargout == 0) || (verbose > 0)
    gsp_plot_graph(G);
    t = 'Spectral layout with';
    for i = 1:length(indices)-1
        t = strcat(t, ' U(', num2str(indices(i)),'),');
    end
    t = strcat(t, ' U(', num2str(indices(end)),')');
    title(t);
end
