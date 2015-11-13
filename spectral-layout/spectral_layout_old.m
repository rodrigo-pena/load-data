function [ G, U ] = spectral_layout_old( G, p, method, TOL, MAX_ITER, verbose )
% SPECTRAL_LAYOUT_OLD computes new coordinates for the nodes of graph G 
% based on an approximation of U(:,2),...,U(:,p), the top (nondegenerate) 
% eigenvectors of D^(-1)*G.W
%
% G is a graph, G.W is its weighted adjacency matrix, and D = sum(G.W, 2) 
% is its diagonal degree matrix. The coordinates of the vertices in G are 
% set to U(:,2:p). If these eigenvectors are already computed, or if you 
% need and exact computation of them, use spectral_layout.m
%
%   Input parameters:
%         G         : Structure containing graph information (see GSPBox
%                     docs). The graph Laplacian eigenvectors should be
%                     accessible at G.U, and the respective eigenvalues at
%                     G.e.
%         p         : Number specifying the index of the topmost
%                     eigenvector to compute.
%         method    : - 'standard': use standard eigenvectors of G.L, the
%                       graph Laplacian.
%                     - 'deg-normalized': use generalized eigenvectors of
%                       (G.L, D), where D is the diagonal degree matrix.
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
%   Example:::
%         G = gsp_bunny(); % Requires GSPBox
%         figure();
%         subplot(121)
%         gsp_plot_graph(G);
%         title('Original Coordinates')
%         subplot(122)
%         G2 = spectral_layout_2(G, 3);
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
% Date: 10 Nov 2015
% Testing: test/test_spectral_layout.m

%% Parse input

% p
p = round(p); % Make sure p is an integer

assert((p <= G.N) && (p > 2), ...
    'p should be an integer in the set the set {3,...,G.N}');

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

% TOL
if (nargin < 4) || (isempty(TOL))
    TOL = 1e-6;
end

% MAX_ITER
if (nargin < 5) || (isempty(MAX_ITER))
    MAX_ITER = 1000;
end

% verbose
if (nargin < 6) || (isempty(verbose))
    verbose = 0;
end

%% Compute the eigenvectors
epsilon = TOL; % Tolerance
U = ones(G.N, p) / sqrt(G.N); % Eigenvector column matrix

% Power-iteration
for k = 2:p
    U_hat = rand(G.N, 1); % Random initialization
    U_hat = U_hat / norm(U_hat);
    
    change = 0;
    niter = 0;
    
    while (change < 1 - epsilon) && (niter <= MAX_ITER)
    
    U(:,k) = U_hat;    
        
    % D-Orthogonalize against previous eigenvectors
    for l = 1:k-1
        U(:,k) = U(:,k) - ...
            ((U(:,k)' * D * U(:,l))/(U(:,l)' * D * U(:,l))) * U(:,l);
    end
    
    % Multiply with (1/2)*(I + D^(-1)*G.W)
    U_hat = (1/2) * (eye(G.N) + D\G.W) * U(:,k);
    
    U_hat = U_hat / norm(U_hat); % Normalization
    
    change = (U_hat' * U(:,k)); 
    niter = niter + 1;
    
    end
    
    U(:,k) = U_hat;
    
end

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
