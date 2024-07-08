function [U_sigma, sigma] = match_eigs(U,V)
%match_eigs matches two sets of orthonormal eigenvectors
% 
%   Command: [U_sigma, sigma] = match_eigs(U,V)
% 
%   U and V are N x N matrices with eigenvectors along the columns.
% 
%   Outputs:
%       U_sigma: N x N matrix of eigenvectors from U matched to V
%       sigma: N x 1 vector corresponding to best permutation.
%       U_sigma == U(:,sigma)


% Ensure U and V are square and of the same size
assert(all(size(U) == size(V)), 'Matrices U and V must be the same size and square.');

% Calculate the cost matrix as 1 minus the absolute value of the dot products
costMatrix = 1 - abs(U' * V);

% Find minimum cost matching using the matchpairs function
sigma = matchpairs(costMatrix, 1000);
sigma = sigma(:,1);

% Arrange U according to sigma to get U_sigma
U_sigma = U(:, sigma);

end