function [ V ] = orthonormalize( U, D )
% ORTHONORMALIZE orthonormalizes a set of vectors. It also orthogonalizes
% the vectors against 1_n (the n-dimensional vector of 1s).
%
%   Input parameters:
%         U         : A n-by-m matrix of column vectors
%         D         : (Optional) A n-by-n matrix. Use when
%                     D-orthonormalization of the columns of U is needed.
%                     Default: eye(n);
%
%   Output parameters:
%         V         : The (D-)orthonormalization of the columns of U.
%
%   Example:
%         n = 200;
%         m = 100;
%         U = rand([n, m]);
%         V = orthonormalize(U);
%
%   See also: hde.m
%
%   Requires:
%
%   References:
%   [1]	Y. Koren, "Drawing Graphs by Eigenvectors," Computers and
%       Mathematics with Applications, vol. 49, pp. 1867-1888, 2005.
%
% Author: Rodrigo Pena
% Date: 13 Nov 2015
% Testing:

%% Parse input

n = size(U,1);
m = size(U,2);

assert(n >= m)

if (nargin < 2) || isempty(D)
    D = eye(n);
end

assert((size(D,1) == n) && (size(D,2) == n), ...
    'D must be a n-by-n matrix, where n = length(U(:,1))');

%% Gram-Schmidt orthonormalization

% Initialization
V = [ones(n,1), U]; % Concatenate 1_n with U
V(:,1) = V(:,1)./norm(V(:,1), 2);
epsilon = 1e-3;

for i = 2:m+1
    % D-orthogonalize
    V(:, i:end) = V(:, i:end) - ...
        repmat((V(:, i:end)' * D * V(:, i - 1))', [n, 1]) .* ...
        repmat(V(:, i - 1), [1, m + 2 - i]);

    if norm(V(:,i),2) < epsilon % A linearly dependent vector
        V(:,i) = zeros([n, 1]);
    else
        % Normalize
        V(:,i) = V(:,i)./norm(V(:,i), 2);
    end
end

V = V(:,2:end); % Remove 1_n

end
