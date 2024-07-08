function output = analyze_network(G, option)
%analyze_network takes in a network and outputs standard graph analyses
%   Usage: output = analyze_network(G,option)
%   options: 'local','global','both'

if ~exist('option','var')
    option = 'both';
end

[n, m] = size(G);

if n ~= m
    error('Input matrix must be square!')
end

G = G.*~eye(n);

if ~issymmetric(G)
%     warning('Making network symmetric!');
    G = triu(G,1) + triu(G,1)';
end

output = struct('global',[],'local',[]);

G(isnan(G)) = 0;
if ismember(0, sum(G,1))
    output.warning = ['Empty ROIs: ', num2str(find(~sign(sum(G,1))))];
end

Dist_G = (1./G);
Dist_G(isnan(Dist_G)) = 0; Dist_G(isinf(Dist_G)) = 0;

% % % Params:
gamma = 1;

[core_structure, coreness_q] = core_periphery_dir(G, gamma);
[charpathlen,~,ecc,radius,diameter] = charpath(Dist_G,0,0);

% [motif_intensity, motif_coherence, motif_frequency] = motif3struct_wei(weight_conversion(G, 'normalize'));


% % % LOCAL MEASURES:
if strcmp(option,'both') || strcmp(option,'local')
    output.local.node_strength = strengths_und_sign(G);
    % output.local.topological_overlap = gtom(G, 1);
    output.local.core_structure = core_structure;
    output.local.node_eccentricity = ecc;
    % output.local.mean_first_pass_time = mean_first_passage_time(G);
    output.local.betweenness_centrality = betweenness_wei(Dist_G)./( (n-1)*(n-2) );
    output.local.eig_centrality = eigenvector_centrality_und(G);
    output.local.pagerank_centrality = pagerank_centrality(G, 0.85);
end


% % % Global MEASURES:
if strcmp(option,'both') || strcmp(option,'global')
    output.global.assortativity = assortativity_wei(G,0);
    output.global.coreness_q = coreness_q;
    output.global.char_path_length = charpathlen;
    output.global.radius = radius;
    output.global.diameter = diameter;
    output.global.small_world = smallworldness(G,'all');
end



if ismember(-1, unique(sign(G(:))))
    G_pos = G;
    G_pos(G<0) = 0;
    
    [M, Q] = community_louvain(G, gamma,[],'negative_sym');
    output.global.modularity_Q = Q;
    output.local.moularity_assignments = M;
    [output.local.clustering_coef_pos, output.local.clustering_coef_neg, output.global.clustering_coef_pos_mean, output.global.clustering_coef_neg_mean] = clustering_coef_wu_sign(G);
    
    % % % LOCAL MEASURES:
    if strcmp(option,'both') || strcmp(option,'local')
        [output.local.local_assortativity_pos, output.local.local_assortativity_neg] = local_assortativity_wu_sign(G);
        output.local.participation_coef = participation_coef_sign(G,M);
        [output.local.diversity_coef_pos, output.local.diversity_coef_neg] = diversity_coef_sign(G,M);
        output.local.local_efficiency = efficiency_wei(weight_conversion(G_pos,'normalize'),2);
    end
    
    % % % GLOBAL MEASURES:
    if strcmp(option,'both') || strcmp(option,'global')
    output.global.efficiency = efficiency_wei(G_pos);
    output.global.transitivity = transitivity_wu(weight_conversion(G_pos,'normalize'));
    
    RC_curve = rich_club_wu(G_pos);
    RC_curve(isnan(RC_curve)) = 0;
    RC_AUC = trapz(RC_curve);
    output.global.rich_club_curve = RC_curve;
    output.global.rich_club_AUC = RC_AUC;
    
    [~,IOTA]=quasi_idempotence(G_pos);
    output.global.quasi_idempotence_IOTA1 = IOTA(1);
    output.global.quasi_idempotence_IOTAinf = IOTA(end);
    end

else

[M, Q] = community_louvain(G, gamma);
output.global.modularity_Q = Q;
output.local.moularity_assignments = M;
output.local.clustering_coef = clustering_coef_wu(G);
output.global.clustering_coef_mean = mean(output.local.clustering_coef);

% % % LOCAL MEASURES:
    if strcmp(option,'both') || strcmp(option,'local')
    output.local.local_efficiency = efficiency_wei(weight_conversion(G,'normalize'),2);
    output.local.local_assortativity = local_assortativity_wu_sign(G);
    output.local.participation_coef = participation_coef(G,M);
    output.local.diversity_coef = diversity_coef_sign(G,M);
    end

    % % % GLOBAL MEASURES:
    if strcmp(option,'both') || strcmp(option,'global')
    output.global.efficiency = efficiency_wei(G);
    output.global.transitivity = transitivity_wu(weight_conversion(G,'normalize'));
    
    RC_curve = rich_club_wu(G);
    RC_curve(isnan(RC_curve)) = 0;
    RC_AUC = trapz(RC_curve);
    output.global.rich_club_curve = RC_curve;
    output.global.rich_club_AUC = RC_AUC;

    [~,IOTA]=quasi_idempotence(G);
    output.global.quasi_idempotence_IOTA1 = IOTA(1);
    output.global.quasi_idempotence_IOTAinf = IOTA(end);
    end

end

end

function out = smallworldness(G,method, nrandomizations)
%smallworldness computes the small worldness of a network G.
%   Usage: out = smallworldness(G,method, nrandomizations)
% 
%   Inputs: 
%       G: a symmetric matrix. 
%           Note that negative weights in G are set to zero.
%       method: options include the following...
%           'sigma': (C/C_r) / (L/L_r) where C is the clustering coef, L is
%           the characteristic path length and _r denotes a
%           degree-preserving randomized network.
%           'omega': 1-abs((L_r/L) - (C/C_l)) where C_l is the clustering coef for
%           an equivalent lattice network.
%           'SWI': ((L - L_l) / (L_r - L_l)) * ((C - C_r) / (C_l - C_r))
%           'Phi': 1 - sqrt(( ((C_l - C)/ (C_l - C_r))^2 + ((L - L_r)/ (L_l - L_r))^2 )/2)
%       nrandomizations: the number of randomized networks to generate for
%       the computation. Values from randomizations are averaged before the
%       final computation.
% 
%   Output is the small-worldness defined by the method.


if ~exist('nrandomizations','var')
    nrandomizations = 1;
end

if ~exist('method','var')
    method = 'sigma';
end

[N, m] = size(G);

if N ~= m
    error('Input matrix must be square!')
end

G = G.*~eye(N);
G(G<0)=0;
G(isnan(G)) = 0;
G = G./mean(nonzeros(triu(G,1)));

if ~issymmetric(G)
%     warning('Making network symmetric!');
    G = triu(G,1) + triu(G,1)';
end

Dist_G = (1./G);
Dist_G(isnan(Dist_G)) = 0; 
Dist_G(isinf(Dist_G)) = 0;

C = mean(clustering_coef_wu(G));
L = charpath(Dist_G,0,0);

rng('shuffle','twister');
C_r = zeros(nrandomizations,1);
L_r = zeros(nrandomizations,1);
for i = 1:nrandomizations
    G_rand = null_model_und_sign(G,2);
    C_r(i) = mean(clustering_coef_wu(G_rand));
    Dist_G_rand = 1./G_rand; 
    Dist_G_rand(isnan(Dist_G_rand))=0;
    Dist_G_rand(isinf(Dist_G_rand))=0;
    L_r(i) = charpath(Dist_G_rand,0,0);
end
C_r = mean(C_r);
L_r = mean(L_r);


switch method
    case 'sigma'
        out = (C/C_r) / (L/L_r);
    case 'omega'
        % % % Define Lattice
        K = sum(triu(weight_conversion(G,'binarize')),'all');
        Lattice_network = makelatticeCIJ(N,K);
        G_nonzeros = nonzeros(triu(G));
        Lattice_network(find(Lattice_network==1)) = G_nonzeros(randi(numel(G_nonzeros), size(find(Lattice_network==1))));
        Lattice_network = triu(Lattice_network,1) + triu(Lattice_network,1)';
        C_l = mean(clustering_coef_wu(Lattice_network));
        out = 1 - abs((L_r/L) - (C/C_l));
    case 'SWI'
        % % % Define Lattice
        K = sum(triu(weight_conversion(G,'binarize')),'all');
        Lattice_network = makelatticeCIJ(N,K);
        G_nonzeros = nonzeros(triu(G));
        Lattice_network(find(Lattice_network==1)) = G_nonzeros(randi(numel(G_nonzeros), size(find(Lattice_network==1))));
        Lattice_network = triu(Lattice_network,1) + triu(Lattice_network,1)';
        C_l = mean(clustering_coef_wu(Lattice_network));
        Dist_Lattice = 1./Lattice_network;
        Dist_Lattice(isnan(Dist_Lattice))=0;
        Dist_Lattice(isinf(Dist_Lattice))=0;
        L_l = charpath(Dist_Lattice);
        out = ((L - L_l) / (L_r - L_l)) * ((C - C_r) / (C_l - C_r));

    case 'phi'
        K = sum(triu(weight_conversion(G,'binarize')),'all');
        Lattice_network = makelatticeCIJ(N,K);
        G_nonzeros = nonzeros(triu(G));
        Lattice_network(find(Lattice_network==1)) = G_nonzeros(randi(numel(G_nonzeros), size(find(Lattice_network==1))));
        Lattice_network = triu(Lattice_network,1) + triu(Lattice_network,1)';
        C_l = mean(clustering_coef_wu(Lattice_network));
        Dist_Lattice = 1./Lattice_network;
        Dist_Lattice(isnan(Dist_Lattice))=0;
        Dist_Lattice(isinf(Dist_Lattice))=0;
        L_l = charpath(Dist_Lattice);
        out = 1 - sqrt(( ((C_l - C)/ (C_l - C_r))^2 + ((L - L_r)/ (L_l - L_r))^2 )/2);

    case 'all'
        % % % Define Lattice
        K = sum(triu(weight_conversion(G,'binarize')),'all');
        Lattice_network = makelatticeCIJ(N,K);
        G_nonzeros = nonzeros(triu(G));
        Lattice_network(find(Lattice_network==1)) = G_nonzeros(randi(numel(G_nonzeros), size(find(Lattice_network==1))));
        Lattice_network = triu(Lattice_network,1) + triu(Lattice_network,1)';
        C_l = mean(clustering_coef_wu(Lattice_network));
        Dist_Lattice = 1./Lattice_network;
        Dist_Lattice(isnan(Dist_Lattice))=0;
        Dist_Lattice(isinf(Dist_Lattice))=0;
        L_l = charpath(Dist_Lattice);

        out.sigma = (C/C_r) / (L/L_r);
        out.omega = 1 - abs((L_r/L) - (C/C_l));
        out.SWI = ((L - L_l) / (L_r - L_l)) * ((C - C_r) / (C_l - C_r));
        out.phi = 1 - sqrt(( ((C_l - C)/ (C_l - C_r))^2 + ((L - L_r)/ (L_l - L_r))^2 )/2);
end

end