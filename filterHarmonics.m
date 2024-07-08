function FilteredNetwork = filterHarmonics(A, weights, sigma, rm_neg)
%filterHarmonics filters a network based on its Laplacian harmonics.
%   Usage: FilteredNetwork = filterHarmonics(A, weights, sigma)
%          FilteredNetwork = filterHarmonics(cat(3,U,diag(ev)), weights, sigma)
% 
%   Inputs:
%       A: a network's adjacency matrix (must be square and symmetric).
%           OR
%          a 3D matrix where A(:,:,1) are Laplacian eigenvectors and 
%          A(:,:,2) is a diagonal matrix of the associated eigenvalues.
% 
%       weights: a vector the same length as A defining one of two possible
%       filtring schemes
%           option 1: categorical weighting. Divides A into N networks
%           based on the interger harmonic assignments in weights.
%           unique(weights) must be equal to 1:N
%           option 2: continuous weighting. Reutrns 1 network where
%           harmonics are weighted based on the weights vector.
%           mean(weights) must be equal to 1.
% 
%       sigma: (optional) a vector defining a permutation for harmonics
%       from ascending order to another ordering.
%       
%       rm_neg: (bool) Removes remaining negative weights in filtered networks.
%           (default: true)
% 
%   Outputs:
%       FilteredNetwork: a network the same number of rows/columns as A,
%       except filtered according to the above. 
%       If weighting was categorical, FilteredNetworks has each category 
%       along the 3rd dimension.
%       If weighting was continuous, only one network is returned.


[sz] = size(A);

if ~isequal(sz(1), sz(2))
    error('Input matrix is not square!');
end

nroi = sz(1);

if length(weights) ~= nroi
    error('The weight vector should have the same number of rows as A.')
end

if ~exist('rm_neg','var') || isempty(rm_neg)
    rm_neg = true;
end

if length(sz) == 3 && sz(3) == 2
    U = A(:,:,1);
    ev = diag(A(:,:,2));

elseif issymmetric(A)

    if nargin < 3 || isempty(sigma)
        sigma = 1:nroi;
    end
    
    D = sum(A, 1);
    nroi = length(D);
    sqrt_D_inv = diag(1./sqrt(D));
    sqrt_D = diag(sqrt(D));
    
    L = eye(nroi) - sqrt_D_inv * A * sqrt_D_inv ;
    
    [U, ev] = eig(L);
    ev = diag(ev);
    [~, ii] = sort(ev, 'ascend');
    ev = ev(ii);
    ev = ev(sigma);
    U = U(:,ii);
    U = U(:, sigma);

else
    error('Input matrix must be either 3D with Eigenvectors and Eigenvalues or a symmetric matrx.')
end

Harmonic_stack = zeros(nroi, nroi, nroi);
% D_stack = zeros(nroi, nroi);

if mean(weights) < 1.01
    
    for i = 1:nroi
        Harmonic_stack(:,:,i) = weights(i) * ev(i) * (U(:,i) * U(:,i)');
        % D_stack(:,i) = diag(Harmonic_stack(:,:,i));
    end

    % D_trunc = diag(sum(D_stack,2));
    % sqrt_D = diag(sqrt(diag(D_trunc)));
    FilteredNetwork = -1*(sqrt_D * sum(Harmonic_stack, 3) * sqrt_D).*~eye(nroi);

else

    n = unique(weights);

    if isequal(n, 1:length(n))
        error('Categorical weights must be ordinal 1:n');
    end

    FilteredNetwork = zeros(nroi, nroi, length(n));

    Harmonic_stack = zeros(nroi, nroi, nroi);
    % D_stack = zeros(nroi, nroi);

    for i = 1:nroi
        Harmonic_stack(:,:,i) = ev(i) * (U(:,i) * U(:,i)');
        % D_stack(:,i) = diag(Harmonic_stack(:,:,i));
    end

    for c = 1:length(n)

        % idx = find(weights == c);

        % D_trunc = diag(sum(D_stack(:,idx),2));
        % D_trunc = diag(D)*diag(sum(D_stack(:,idx),2));

        % sqrt_D = diag(sqrt(diag(D_trunc)));
        

        FilteredNetwork(:,:,c) = -1*(sqrt_D * sum(Harmonic_stack(:,:,weights == c), 3) * sqrt_D).*~eye(nroi);

    end

end

if rm_neg
    FilteredNetwork(FilteredNetwork < 0) = 0;
end


end