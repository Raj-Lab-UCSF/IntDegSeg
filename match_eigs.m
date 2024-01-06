function [U_matched, matched_order, UU_dot] = match_eigs(U1,U2)
%match_eigs will find a best match bewteen two sets of
%eigenvectors, where the first (U1) is truncated to the second's (U2) size. 
% 
%   Command: [U1_matched, , matched_order, UU_dot] = match_eigs(U1, U2)
%   Inputs:
%       U1: N x M matrix of truncated (via rows) eigenvectors to match U2
%       U2: N x N matrix of eigenvectors to match to U1
% 
%   Outputs:
%       U1_matched: N x N matrix of eigenvectors that pseudo-optimally
%       matched U2
%       matched_order: N x 1 vector corresponding to the ordering of
%       eigenvectors in U1. U_matched == U1(:,matched_order)
%       UU_dot: NxM matrix = abs(U1' * U2)

neigs_cut = length(U2);
neigs_full = size(U1,2);

UU_dot = abs(U1' * U2); %projection of natural eigs to non-subcortex eigs

U_matched = zeros(neigs_cut);
matched_order = zeros(neigs_cut,1);

missing_eigs = 1:neigs_cut;
available_eigs = 1:neigs_full;

for i = 1:neigs_cut

    tmp = UU_dot(available_eigs,missing_eigs);
    [~, idx] = max(tmp,[],'all');
    [idx_i, idx_j] = ind2sub(size(tmp),idx);
    U_matched(:,missing_eigs(idx_j)) = U1(:,available_eigs(idx_i));

    matched_order(missing_eigs(idx_j)) = available_eigs(idx_i);
    available_eigs(idx_i) = [];
    missing_eigs(idx_j) = [];

end


end