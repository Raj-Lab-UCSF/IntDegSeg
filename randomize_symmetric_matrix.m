function B = randomize_symmetric_matrix(A)
% RANDOMIZE_SYMMETRIC_MATRIX shuffles the upper triangular part of a
% symmetric matrix and returns a new symmetric matrix with the shuffled
% values.
%
% Input:
%   - A: a symmetric matrix
%
% Output:
%   - B: a new symmetric matrix with the upper triangular part of A
%   randomly shuffled

% Get the size of the matrix
n = size(A,1);

% Extract the upper triangular part of A as a vector
A_vec = A(triu(true(n), 1));

% Shuffle the vector randomly
A_vec_shuffled = A_vec(randperm(length(A_vec)));

% Create a new matrix with the shuffled upper triangular part
B = zeros(size(A));
B(triu(true(n), 1)) = A_vec_shuffled;

% Copy the lower triangular part of B to the upper triangular part
B = triu(B) + triu(B, 1)';

end
