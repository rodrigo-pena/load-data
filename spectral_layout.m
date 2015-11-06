function [ G ] = spectral_layout( G, indices, method, verbose )
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
%
%   Output parameters:
%         G         : Structure containing graph information (see GSPBox
%                     docs), with the nodes' coordinates updated.
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
%   See also:
%
%   Requires: GSPBox (https://lts2.epfl.ch/gsp/)
%
%   References:
%   [1]	Y. Koren, "Drawing Graphs by Eigenvectors," Computers and
%       Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.
%
% Author: Rodrigo Pena
% Date: 6 Nov 2015
% Testing: test/test_spectral_layout.m

%% Parse input

assert(sum(indices>G.N)==0 && sum(indices<1)==0 && length(indices)<=G.N,...
    'indices should not exceed G.N in length');

assert(sum(indices>G.N)==0 && sum(indices<1)==0,...
    'indices should not have elements out of the set {1,...,G.N}');

indices = round(indices); % Make sure the elements in indices are integers

if (nargin < 3) || (isempty(method))
    method = 'standard';
end

assert(strcmp(method, 'standard') || strcmp(method, 'deg-normalized'), ...
    'method should be `standard`, or `deg-normalized`');

if (nargin < 4) || (isempty(verbose))
    verbose = 0;
end

%% Compute new coordinates
switch method
    case 'standard'
        
        if ~isfield(G, 'U') || ~isfield(G, 'e')
            G.U = [];
            G.e = [];
        end
        
        if length(G.e) < max(indices)
            warning('Pre-compute the eigendecomposition of L for speed.');
            
            if ~issparse(G.L) || (G.N <= 3000)
                try
                    [G.U, G.e] = eig(full(G.L));
                catch
                    [G.U, G.e, ~] = svd(full(G.L));
                end
            else
                try
                    [G.U, G.e] = eigs(G.L, [], max(indices), 'SM');
                catch
                    [G.U, G.e, ~] = svds(G.L, max(indices), 0);
                end
            end
            
            G.e = diag(G.e);
            [G.e, idx] = sort(G.e);
            G.U = G.U(:,idx);
            
            % Choose non-negative representation
            signs = sign(G.U(1,:));
            signs(signs==0)=1;
            G.U = G.U * diag(signs);
            
        end
                
        G.coords = G.U(:, indices);
        G = gsp_graph_default_plotting_parameters(G);
        minVec = min(G.coords, [], 1);
        maxVec = max(G.coords, [], 1);
        G.plotting.limits = reshape([minVec; maxVec], 1, []);
        
    case 'deg-normalized'
        
        if ~isfield(G, 'U_gen') || ~isfield(G, 'e_gen')
            G.U_gen = [];
            G.e_gen = [];
        end
        
        if length(G.e_gen) < max(indices)
            
            warning('Pre-compute the generalized eigendecomposition of (L,D) for speed.');
            
            if ~isfield(G, 'd')
                G.d = sum(G.W, 2);
            end
            D = diag(G.d);
            
            if ~issparse(G.L) || (G.N <= 3000)
                try
                    [G.U_gen, G.e_gen] = eig(full(G.L),full(D));
                catch
                    [G.U_gen, G.e_gen, ~] = svd(full(D\G.L));
                end
            else
                try
                    [G.U_gen, G.e_gen] = eigs(G.L, D, max(indices), 'SM');
                catch
                    [G.U_gen, G.e_gen, ~] = svds(D\G.L, max(indices), 0);
                end
            end
            
            G.e_gen = diag(G.e_gen);
            [G.e_gen, idx] = sort(G.e_gen);
            G.U_gen = G.U_gen(:,idx);
            
            % Choose non-negative representation
            signs = sign(G.U_gen(1,:));
            signs(signs==0)=1;
            G.U_gen = G.U_gen * diag(signs);
            
        end
        
        G.coords = G.U_gen(:, indices);
        G = gsp_graph_default_plotting_parameters(G);
        minVec = min(G.coords, [], 1);
        maxVec = max(G.coords, [], 1);
        G.plotting.limits = reshape([minVec; maxVec], 1, []);
end

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
