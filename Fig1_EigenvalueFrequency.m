%% Figure 1: Eigenvalue Frequency Script

%% Setup:

scale = 200;

filename = sprintf('MICA_schaefer%d_SCFC_struct.mat',scale);
% filename = sprintf('MICA_HCP_schaefer%d_SCFC_struct.mat',scale);
load(filename);

n_subjects = length(MICA);
nroi = length(MICA(1).SC);

SC_all = zeros(nroi,nroi, n_subjects);
SC_ev_all = zeros(n_subjects,nroi);

for i = 1:n_subjects
    SC = MICA(i).SC;
    SC = SC./norm(SC,'fro');
    SC_all(:,:,i) = SC;
    [~, ~, SC_ev_all(i,:)] = graph_laplacian(SC, 'normalized');
end
% clear MICA;

SC_consensus = mean(SC_all,3);
[SC_L, SC_U_consensus, SC_ev_consensus] = graph_laplacian(SC_consensus, 'normalized');


sc_color = [1 0.25 0.5];
roughness_color = [1.0000, 0.7969, 0.3984]; %gold
zeroX_color =  [0.933, 0.475, 0.176]; %'bright orange'
support_color = [0.9219, 0.4375, 0.3867]; %coral
mincut_color = [0.9531, 0.6875, 0.5156]; %apricot  

%% Eigenvalue Plot:
FIG = figure();
hold on;
p = plot(1:length(SC_ev_consensus), SC_ev_consensus);
p.Color = sc_color;
p.LineWidth = 5;
xlim([0, scale+20]);
ylim([0, 1.6]);

p1 = plot_iqr([], SC_ev_all, 'median', [0, 0, 0], true);
p1.LineWidth = 2;

title('Eigenvalue Plot consensus versus Subjects');
xlabel('Eigenmode Index');
ylabel('Eigenvalues');
box on;
set(FIG, 'Position',[1,49,1920,955]);

%% Visualize Harmonics (Requires BrainNet Viewer)

idx = 4;
surface_node_idx = ones(200,1);
config_file = 'BrainNet_Config_full_surface_and_node.mat';
save_name = ['Consensus_Harmonic_' num2str(idx)];

BNV_surface_and_nodes(SC_U_consensus(15:end, idx), 'Schaefer_200.nii.gz', [], surface_node_idx, config_file, 'Ch2',save_name);

%% % % Find Zerocrossrate, entropy, min-cut-max-flow, and support

load(sprintf('Schaefer%d_Adjacency.mat',scale));

Adja_norm = Adja./norm(Adja,'fro');

LH = (14+1):(14+(scale/2));
RH = (14+(scale/2)+1):(scale+14);

Adja_L = Adja(LH, LH);
Adja_R = Adja(RH, RH);
Adja_L_norm = Adja_L./norm(Adja_L,'fro');
Adja_R_norm = Adja_R./norm(Adja_R,'fro');
L_adja_L = graph_laplacian(Adja_L_norm,'normalized');
L_adja_R = graph_laplacian(Adja_R_norm,'normalized');
SC_ev_roughness = mean([vecnorm(L_adja_L*SC_U_consensus(LH,:))', vecnorm(L_adja_R*SC_U_consensus(RH,:))'],2);

% threshold:
thresh = 1e-3;

SC_ev_zeroX = zeros(nroi,1);
SC_ev_sparsity = zeros(nroi,1);
SC_ev_mincut = ev_zeroXings(SC_consensus, SC_U_consensus, thresh);

for i = 1:nroi

    evec = SC_U_consensus(:,i);

    evec_thresh = evec;
    evec_thresh(abs(evec_thresh) <=thresh) = 0;

    SC_ev_zeroX(i) = zerocrossrate(evec_thresh);
    
    evec_thresh(abs(evec_thresh) >thresh) = 1;

    SC_ev_sparsity(i) = (nroi - sum(evec_thresh))/nroi;

end

%% Pannel B: Scatter Plots "Ordered" Frequency

sz = 200;

FIG = figure();
scatter(SC_ev_consensus, SC_ev_sparsity,sz,'filled', 'MarkerFaceColor', support_color, 'MarkerEdgeColor', 'black');
box on;
xlim([0, 1.5]);
pbaspect([3.63, 1.81, 1]);
set(FIG, 'Position',[1,49,1920,955]);

FIG = figure();
scatter(SC_ev_consensus, SC_ev_zeroX,sz,'filled', 'MarkerFaceColor', zeroX_color, 'MarkerEdgeColor', 'black');
box on;
xlim([0, 1.5]);
ylim([0,0.7]);
pbaspect([3.63, 1.81, 1]);
set(FIG, 'Position',[1,49,1920,955]);

FIG = figure();
scatter(SC_ev_consensus, SC_ev_mincut,sz,'filled', 'MarkerFaceColor', mincut_color, 'MarkerEdgeColor', 'black');
box on;
xlim([0, 1.5]);
pbaspect([3.63, 1.81, 1]);
set(FIG, 'Position',[1,49,1920,955]);

FIG = figure();
scatter(SC_ev_consensus, SC_ev_roughness,sz,'filled', 'MarkerFaceColor', roughness_color, 'MarkerEdgeColor', 'black');
box on;
xlim([0, 1.5]);
pbaspect([3.63, 1.81, 1]);
set(FIG, 'Position',[1,49,1920,955]);

