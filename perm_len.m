function [element_len, len_v] = perm_len(v)
%perm_len outputs a vector that is the element-wise permutation length of
%the input, v, assuming the input is a permutation of 1:length(v).
%   [element_len, len_v] = perm_len(v)
%   Input:
%       v: an Nx1 vector describing the permutation 1:N -> v
%   Output:
%       element_len: an Nx1 vector where each element is its length in the
%       permutation
%       len_v == sum(element_len)/2: the total length of the permutation.

N = length(v);
element_len = zeros(N,1);

for i = 1:N

    if v(i) == i
        continue
    else
        element_len(i) = abs(v(i) - i);
    end

end

len_v = sum(element_len)/2;

end