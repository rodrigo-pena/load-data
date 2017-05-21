function [ G, U ] = spectral_layout_approx(G, p, m, method, TOL, MAX_ITER, ...
    verbose)
% SPECTRAL_LAYOUT_APPROX computes new coordinates for the nodes of graph G 
% based on a subspace-constrained approximation of U(:,2),...,U(:,p), the 
% top (nondegenerate) eigenvectors of B = mu*I - X' * D\G.L * X.
%
% G is a graph, and G.L is its Laplacian. The coordinates of the vertices
% in G are set to U(:,2:p). The columns of the matrix X span the subspace
% to which the coordinates of the nodes of G are constrained.
%
%   Input parameters:
%         G         : Structure containing graph information (see GSPBox
%                     docs). The graph Laplacian eigenvectors should be
%                     accessible at G.U, and the respective eigenvalues at
%                     G.e.
%         p         : Number specifying the index of the topmost
%                     eigenvector to compute.
%         m         : Dimension of the column space of X. Default:
%                     min(G.N, 50).
%         method    : A string. Default: 'standard'
%                       - 'standard': Standard opt. problem (cf. [1])
%                       - 'deg-normalized': Degree-normalized opt. problem
%                         (cf. [1])
%         TOL       : Tolerance of the approximation. 1 - TOL is the cosine
%                     of the angle between the estimated eigenvector at
%                     each iteration. Default: 1e-6.
%         MAX_ITER  : Max number of power-iterations for the approximation
%                     of each eigenvector. Default: 1000.
%         verbose   : Level of information display. Default: 0.
%
%   Output parameters:
%         G         : Structure containing graph information (see GSPBox
%                     docs), with the nodes' coordinates updated.
%         U         : The computed eigenvectors.
%
%   Example:
%         G = gsp_bunny(); % Requires GSPBox
%         figure();
%         subplot(121)
%         gsp_plot_graph(G);
%         title('Original Coordinates')
%         G2 = spectral_layout_approx(G, 4, 100);
%         subplot(122)
%         gsp_plot_graph(G2);
%         title('Spectral Coordinates')
%
%   See also: spectral_layout.m
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   References:
%   [1]	Y. Koren, "Drawing Graphs by Eigenvectors," Computers and
%       Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.
%
% Author: Rodrigo Pena
% Date: 13 Nov 2015
% Testing: test/test_spectral_layout.m

%% Parse input

% G
assert(isfield(G, 'W'));
if ~isfield(G, 'L'); G.L = sum(G.W,2) - G.W; end

% p
p = round(p); % Make sure p is an integer
assert((p <= G.N) && (p > 2), ...
    'p should be an integer in the set the set {3,...,G.N}');

% m
if (nargin < 3)
    m = 50;
end

% method
if (nargin < 4) || (isempty(method))
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

% TOL
if (nargin < 5) || (isempty(TOL))
    TOL = 1e-6;
end

% MAX_ITER
if (nargin < 6) || (isempty(MAX_ITER))
    MAX_ITER = 1000;
end

% verbose
if (nargin < 7) || (isempty(verbose))
    verbose = 0;
end

%% Compute the subspace X

X = hde(G, m);
X = orthonormalize(X, D);

XLX = (X' * ((D\G.L) * X));
dXLX = diag(XLX);
mu = max(dXLX + sum((XLX - diag(dXLX)), 2)); % Gershgorin bound

%% Compute the top eigenvectors of B = mu*I - X' * D\G.L * X
epsilon = TOL; % Tolerance
I = eye(m);
B = mu * I - XLX;
XDX = X' * D * X;
U = ones(m, p) / sqrt(m); % Eigenvector column matrix


% Power-iteration
for k = 2:p
    U_hat = rand(m, 1); % Random initialization
    U_hat = U_hat / norm(U_hat);
    
    change = 0;
    niter = 0;
    
    while (change < 1 - epsilon) && (niter <= MAX_ITER)
        
        U(:,k) = U_hat;
        
        % D-Orthogonalize against previous eigenvectors
        for l = 1:k-1
            U(:,k) = U(:,k) - ...
                ((U(:,k)' * XDX * U(:,l))/(U(:,l)' * XDX * U(:,l))) * U(:,l);
        end
        
        % Multiply by (1/2)*(I + B)
        U_hat = (1/2) * (I + B) * U(:,k);
        
        U_hat = U_hat ./ norm(U_hat); % Normalization
        
        change = (U_hat' * U(:,k));
        niter = niter + 1;
        
    end
    
    U(:,k) = U_hat;
    
end

%% Revert the eigenvectors to the original space frame
U = X*U;

%% Assign the new coordinates
G.coords = U(:,2:p);
G = gsp_graph_default_plotting_parameters(G);
minVec = min(G.coords, [], 1);
maxVec = max(G.coords, [], 1);
G.plotting.limits = reshape([minVec; maxVec], 1, []);

%% Output result
if (nargout == 0) || (verbose > 0)
    gsp_plot_graph(G);
    t = 'Spectral layout with';
    for i = 2:p-1
        t = strcat(t, ' U(', num2str(i),'),');
    end
    t = strcat(t, ' U(', num2str(p),')');
    title(t);
end
