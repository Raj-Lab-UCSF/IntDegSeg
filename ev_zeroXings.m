function [zeroX] = ev_zeroXings(SC, U, thresh)
%ev_zeroXings Computes the zero crossings for all eigenvectors in SC
%   This computes the min-cut/max-flow of the graph. 
%   Usage: min_cut = ev_zeroXings(SC, U, thresh)
%       

    nroi = length(SC);

    if ~exist('thresh','var')
        thresh = zeros(1,nroi);
    elseif thresh >= 1
        thresh = prctile(abs(U).^2,thresh);
    else
        thresh = ones(1,nroi).*thresh;
    end

    zeroX = zeros(nroi,1);

    combos = nchoosek(1:nroi,2);

    for k = 1:nroi
        zeroXk = zeros(length(combos),1);

        V = U(:,k);
        V(abs(V).^2 < thresh(k)) = 0;
        for n = 1:length(combos)
    
            i = combos(n,1);
            j = combos(n,2);
    
% % % % From Huang 2016 article:
            sgn_vi_x_vj = sign(V(i) * V(j));

            if sgn_vi_x_vj == -1
                zeroXk(n) = SC(i,j);
            end

            
        end

        zeroX(k) = 0.5*sum(zeroXk);
    end

end