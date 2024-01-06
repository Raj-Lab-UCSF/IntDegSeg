%% Figure 2: Inter-Subject Agreement Analsysis

%% Inter-Subject Agreement Analsysis

U_Scales_IndVar = struct('scale',[],'UU_diag_matched2tmp',[],'kendall_tau_corr',[], 'all_combo_diagonals',[],'all_combo_kendall_tau_corr',[], 'all_combo_diagonals_matched2TMP',[], 'all_combo_diagonals_matched2IND',[],'permutation_length',[],'L',[], 'matched2IND_permutation_length',[],'matched2IND_L',[]);

scale = 200;

filename = sprintf('MICA_schaefer%d_SCFC_struct.mat',scale);
load(filename);

n_subjects = length(MICA);
nroi = length(MICA(1).SC);

SC_all = zeros(nroi,nroi, n_subjects);
for i = 1:n_subjects
    SC = MICA(i).SC;
    SC = SC./norm(SC,'fro');
    SC_all(:,:,i) = SC;
end
clear MICA;

SC_consensus = mean(SC_all,3);
[~, SC_U_consensus, SC_ev_consensus] = graph_laplacian(SC_consensus, 'normalized');


U_all = zeros(nroi,nroi,n_subjects);
U_all_matched = zeros(nroi,nroi,n_subjects);
matched_order_all = zeros(nroi,n_subjects);
ev_all_matched = zeros(nroi, n_subjects);
UU_diag_matched2tmp = zeros(n_subjects,nroi);
kendall_tau_corr = zeros(n_subjects,1);
permutation_length = zeros(nroi,n_subjects);
L = zeros(n_subjects,1);

for i = 1:n_subjects

    SC = SC_all(:,:,i);
    [~, U_all(:,:,i), ev_sub] = graph_laplacian(SC, 'normalized');

    [U_all_matched(:,:,i), matched_order_all(:,i)] = match_eigs(U_all(:,:,i), SC_U_consensus);

    ev_all_matched(:,i) = ev_sub(matched_order_all(:,i));

    UU_diag_matched2tmp(i,:) = abs(diag(U_all_matched(:,:,i)' * SC_U_consensus));

    kendall_tau_corr(i) = corr(matched_order_all(:,i), [1:nroi]','Type','Kendall');

    [permutation_length(:,i), L(i)] = perm_len(matched_order_all(:,i));
end

combos = nchoosek(1:n_subjects,2);

all_combo_diagonals = zeros(length(combos), nroi);
all_combo_diagonals_matched2TMP = zeros(length(combos), nroi);
all_combo_diagonals_matched2IND = zeros(length(combos), nroi);
all_combo_kendall_tau_corr = zeros(length(combos),1);
matched2IND_permutation_length = zeros(nroi, length(combos));
matched2IND_L = zeros(length(combos),1);

for i = 1:length(combos)

    U1 = U_all(:,:,combos(i,1));
    U2 = U_all(:,:,combos(i,2));

    
% % % % % % No Matching        
    all_combo_diagonals(i,:) = abs(diag(U1' * U2));
%     
% % % % % % Matched Harmonics   
    U1 = U_all_matched(:,:,combos(i,1));
    U2 = U_all_matched(:,:,combos(i,2));
    all_combo_diagonals_matched2TMP(i,:) = abs(diag(U1' * U2));

% % % % % % Best Possible 1-to-1 matching
    [U1_matched2U2, matched2IND_order_all] = match_eigs(U1, U2);
    all_combo_diagonals_matched2IND(i,:) = abs(diag(U1_matched2U2' * U2));
%     
    [matched2IND_permutation_length(:,i), matched2IND_L(i)] = perm_len(matched2IND_order_all);

    all_combo_kendall_tau_corr(i) = corr(matched_order_all(:,combos(i,1)), matched_order_all(:,combos(i,2)),'Type','Kendall');

end

U_Scales_IndVar.scale = scale;
U_Scales_IndVar.UU_diag_matched2tmp = UU_diag_matched2tmp;
U_Scales_IndVar.kendall_tau_corr = kendall_tau_corr;
U_Scales_IndVar.all_combo_diagonals = all_combo_diagonals;
U_Scales_IndVar.all_combo_diagonals_matched2TMP = all_combo_diagonals_matched2TMP;
U_Scales_IndVar.all_combo_diagonals_matched2IND = all_combo_diagonals_matched2IND;
U_Scales_IndVar.matched2IND_permutation_length = matched2IND_permutation_length;
U_Scales_IndVar.matched2IND_L = matched2IND_L;
U_Scales_IndVar.all_combo_kendall_tau_corr = all_combo_kendall_tau_corr;
U_Scales_IndVar.permutation_length = permutation_length;
U_Scales_IndVar.L = L;


%% View Matching
sub = 1;

FIG = figure();
subplot(1,2,1);
imagesc(abs(U_all(:,:,sub)'*SC_U_consensus));
title('Harmonics Before Matching');
colorbar;
axis square;

subplot(1,2,2);
imagesc(abs(U_all_matched(:,:,sub)' * SC_U_consensus));
title('Harmonics Matched to Consensus');
colorbar;
axis square;


%% Figure MICA: Compare subject reliability to each other

dir_figure = 'C://Users//flame//Box//Projects//Eigenmode_Variability//Figures//Experiment1//';



all_combo_diagonals= U_Scales_IndVar.all_combo_diagonals;
all_combo_diagonals_matched2TMP = U_Scales_IndVar.all_combo_diagonals_matched2TMP;
all_combo_diagonals_matched2IND = U_Scales_IndVar.all_combo_diagonals_matched2IND;

nroi = size(all_combo_diagonals,2);

X = 1:nroi;

FIG = figure();
subplot(3,1,1); 
plot_iqr(X, all_combo_diagonals, 'median', [1 0 0.25], true, 0.6);
title('How similar are Individuals to each other? (Unmatched To consensus)');
xlabel('Eigenvector Index');
ylabel('UU Diagonal');
xlim([0, scale+20]);

subplot(3,1,2);
plot_iqr(X, all_combo_diagonals_matched2TMP, 'median', [1 0 0.25], true, 0.6);
title('How similar are Individuals to each other? (Matched To consensus)');
xlabel('Eigenvector Index');
ylabel('UU Diagonal');
xlim([0, scale+20]);

subplot(3,1,3);
plot_iqr(X, all_combo_diagonals_matched2IND, 'median', [1 0 0.25], true, 0.6);
title('How similar are Individuals to each other? (Subs Matched 1-to-1)');
xlabel('Eigenvector Index');
ylabel('UU Diagonal');
xlim([0, scale+20]);

sgtitle(sprintf('Scale: %d',scale));
set(FIG, 'Position', [465,50,937,946]);


%% Permutation Lengths and Kendall Tau:

FIG = figure();
subplot(2,1,1);
plot_iqr(1:(scale+14), permutation_length, 'median', [0 0.4 0.8], true, 0.6);
xlim([0, scale+20]);
title('Permutation Length Per Harmonic');
xlabel('Eigenvector Index');
ylabel('Perm Length');

subplot(2,1,2);
histogram(kendall_tau_corr);
title('Kendalls Tau Permutation Distance');
xlabel('Permutation Distance');
ylabel('Number of Subjects');
set(FIG, 'Position', [465,50,937,946]);

