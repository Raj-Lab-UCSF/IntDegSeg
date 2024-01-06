%% Fig5 RSN and Neurosynth SC Analysis:

% Necessary for Ternary Plots:
addpath ternary_plots/data_plots/
addpath ternary_plots/axes_creation/
addpath ternary_plots/utilities/

%% RSN Analysis:

scale = 200;
filename = sprintf('MICA_schaefer%d_SCFC_struct.mat',scale);
load(filename);
NetworkLabels = readmatrix(sprintf('Schaefer_RSN_Labels/SchaeferSystemLabels_%d.csv',scale));


% Extract SC and make template
n_subjects = length(MICA);
nroi = length(MICA(1).SC);
nt = size(MICA(50).fMRI_TS,1);

SC_all = zeros(nroi,nroi, n_subjects);
SC_U_all = zeros(nroi,nroi, n_subjects);
SC_ev_all = zeros(nroi, n_subjects);

for i = 1:n_subjects
    SC = MICA(i).SC;
    SC = SC./norm(SC,'fro');
    SC_all(:,:,i) = SC;
    [~, SC_U_all(:,:,i), SC_ev_all(:,i)] = graph_laplacian(SC, 'normalized');

end
clear MICA;


SC_template = mean(SC_all,3);
[~, SC_U_template, SC_ev_template] = graph_laplacian(SC_template, 'normalized');


NetworkParticipation_SC_template = NetworkLabels' * abs(SC_U_template).^2; %formal Wavefunction probability distribution function

NetworkParticipation_SC_all = zeros(size(NetworkLabels,2), nroi, n_subjects);


for i = 1:n_subjects

    SC = SC_all(:,:,i);

    [~, SC_U_sub, ~] = graph_laplacian(SC, 'normalized');

    SC_U_sub_matched = match_eigs(SC_U_sub, SC_U_template);
    SC_U_sub_matched_born = (abs(SC_U_sub_matched).^2);
    SC_U_sub_matched_born = SC_U_sub_matched_born./sum(SC_U_sub_matched_born,1); %normalize harmonics for dot product

    NtwrkPart_SC_tmp = NetworkLabels' * SC_U_sub_matched_born;

    NetworkParticipation_SC_all(:,:,i) = NtwrkPart_SC_tmp;

end


%% RSN Group Averaged Heatmap:
Networks = {'Subcortex','Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};

% load('PreComputed_Data\U_Scales_Networks.mat');
% s = 2;
% NetworkParticipation_SC_all = U_Scales_Networks(s).NetworkParticipation_SC_all;

FIG = figure();
NP = NetworkParticipation_SC_all;
NP_mean = mean(NP,3);
% imagesc(NP_mean);

NP_scaled = NP_mean ./ max(NP_mean,[],2);
imagesc(NP_scaled);

colormap hot;
colorbar; clim([0, 1]);
yticklabels(Networks);
title('Yeo Network Heatmap for SC (Mean Across Subjects)');
set(FIG, 'Position',[192,50,1314,946]);

%% RSN Ternary Plot:

% % From EV Regime Consensus:
int_range = 1:26;
deg_range = 27:187;
seg_range = 188:214;


Networks = {'Subcortex','Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
n_subjects = 50;
nnets = length(Networks);
s=2;

% NetworkParticipation_SC_all = U_Scales_Networks(s).NetworkParticipation_SC_all;


pcthresh = 80;
thresh = prctile(NetworkParticipation_SC_all(:),pcthresh);
NetworkParticipation_SC_all_thresh = NetworkParticipation_SC_all;
NetworkParticipation_SC_all_thresh(NetworkParticipation_SC_all < thresh) = 0;

rsn_intdegseg_all = zeros(nnets, 3, 50);

for i = 1:nnets
    for j = 1:n_subjects

        rsn_intdegseg_all(i,1,j) = mean(NetworkParticipation_SC_all_thresh(i,int_range,j));
        rsn_intdegseg_all(i,2,j) = mean(NetworkParticipation_SC_all_thresh(i,deg_range,j));
        rsn_intdegseg_all(i,3,j) = mean(NetworkParticipation_SC_all_thresh(i,seg_range,j));
    end
end

rsn_intdegseg = mean(rsn_intdegseg_all,3);

coords = rsn_intdegseg./sum(rsn_intdegseg,1);
coords = coords./sum(coords,2);
coords = coords(:,[1 3 2]);

% % % YMC Colormap
ycm_colors = zeros(size(coords));  % Initialize colors
ycm_colors(:, 1) = (1 - coords(:, 1)); % Cyan
ycm_colors(:, 2) = (1 - coords(:, 2)); % Yellow
ycm_colors(:, 3) = (1 - coords(:, 3)); % Magenta


coords = rsn_intdegseg./sum(rsn_intdegseg,2);

wlimits = ternary_axes_limits;

FIG = figure();
vgen  = { 'wlimits',       wlimits ,... % Axes will match wlimits ranges
'titlelabels', {'Int','Seg','Deg'}, ... % custom labels
'titlerotation', [0,0,0], ... % Set all titles to horizontal
'link_color', {'tick','title'},... % Link all axes colors
'titleshift',[ -0.18, 0, 0.18; 0.085, -0.11, 0.085 ]... % shift titles
'tickshift', [-0.02, -0.02, 0.0; -0.01,-0.00,0]
};
% Ternary Axes Outline   - Passed directly to plot3()
vout  = { 'LineWidth', 3, 'LineStyle', '-','Color','k'};
% Ternary gridlines  - Passed directly to plot3()
vgrid = { 'LineStyle','--','LineWidth',0.5, 'Color',[0 0 0 0.2] };
% Ternary Tick Line - Passed directly to plot3()
vtick_line = { 'LineStyle', '-', 'LineWidth', 1.5, 'Color',[0 0 0 1.0] };
% Ternary Tick Labels - Passed directly to text()
vtick_label = { 'FontWeight','Bold', 'FontSize', 8 };
% Ternary Axes Label  - Passed directly to text()
vlab  = { 'FontWeight','normal', 'FontSize', 14 };
h = ternary_axes( vgen, vout, vgrid, vtick_line, vtick_label, vlab );
[pH, ~, flat_coords] = ternary_scatter3( wlimits, 'l', coords(:,1), 'b', coords(:,3), coords(:,2));
pH.MarkerEdgeColor = [0, 0, 0];
pH.CData = ycm_colors;
pH.SizeData = 200;
colorbar off;

hold on;
for i = 1:length(Networks)
    text(flat_coords(i,1)+0.015, flat_coords(i,2), flat_coords(i,3), Networks(i))
end

set(FIG, 'Position',[1,49,1920,955]);


%% Neurosynth Meta Analysis Load:
load('Domain_Permutation.mat'); %from excel sheet, saved as .mat for convenience
domain_groupings = domain_labels(domain_perm);

n_subjects = 50;
scale=200;

EigMeta200_all = zeros(146, scale+14, n_subjects);

for sub = 1:n_subjects
    EigMeta200_all(:,:,sub) = readmatrix(sprintf('Neurosynth_Results/Sub%d_Scale%d_absUsquared_LDA400_Metanalysis_Results.csv',sub,scale), 'Range',[2 2]);
end

% % % Analysis without Thresholding the LDA groups:
ngroups = 25;

EigMeta200_all_perm = EigMeta200_all(domain_perm,:,:);
EigMeta200_all_grouped = zeros(ngroups, size(EigMeta200_all_perm,2), n_subjects);
for i = 1:ngroups
    EigMeta200_all_grouped(i,:,:) = mean(EigMeta200_all_perm(find(domain_groupings==i),:,:),1);
end

% Remove Pain from List:
EigMeta200_all_grouped(find(domain_names=='Pain'),:,:) = [];
domain_names(find(domain_names=='Pain')) = [];
ngroups = length(domain_names);

EigMeta200_all_grouped_SC = EigMeta200_all_grouped;



%% Meta Analysis Mean Across Subjects


EigMeta200_all_mean_SC = mean(EigMeta200_all_grouped_SC,3);
EigMeta200_all_mean_SC(EigMeta200_all_mean_SC<0) = 0; %remove negative correlations

FIG = figure();
imagesc(EigMeta200_all_mean_SC); colormap hot; colorbar;

set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);

title('Mean Across Subjects Meta Analysis (SC)');
yticks(1:24); yticklabels(domain_names);
set(FIG, 'Position',[192,50,1314,946]);





%% Ternary Plot Neurosynth

% % From EV Regime Consensus:
int_range = 1:26;
deg_range = 27:187;
seg_range = 188:214;


coord_cmap = jet(24);

domain_names_abbv = strings(24,1);
domain_names_abbv(1) = 'VisPer';
domain_names_abbv(2) = 'MSP';
domain_names_abbv(3) = 'Mot';
domain_names_abbv(4) = 'NumCog';
domain_names_abbv(5) = 'Aud';
domain_names_abbv(6) = 'VisAtn';
domain_names_abbv(7) = 'Awr';
domain_names_abbv(8) = 'Act';
domain_names_abbv(9) = 'VisSpa';
domain_names_abbv(10) = 'CueAttn';
domain_names_abbv(11) = 'WMem';
domain_names_abbv(12) = 'Inh';
domain_names_abbv(13) = 'Lang';
domain_names_abbv(14) = 'CogCon';
domain_names_abbv(15) = 'VerbSem';
domain_names_abbv(16) = 'DMem';
domain_names_abbv(17) = 'ABMem';
domain_names_abbv(18) = 'Pred';
domain_names_abbv(19) = 'SocCog';
domain_names_abbv(20) = 'VisSem';
domain_names_abbv(21) = 'FA';
domain_names_abbv(22) = 'Em';
domain_names_abbv(23) = 'DecM';
domain_names_abbv(24) = 'ImCre';

left_paren = strings(24,1);
left_paren(:) = ' (';
right_paren = strings(24,1);
right_paren(:) = ')';

legend_labels = strcat(domain_names, left_paren, domain_names_abbv, right_paren);

EigMeta200_all_grouped_SC(EigMeta200_all_grouped_SC<0) = 0;

pcthresh = 80;
thresh = prctile(EigMeta200_all_grouped_SC(:),pcthresh);

EigMeta200_all_grouped_SC_thresh = EigMeta200_all_grouped_SC;
EigMeta200_all_grouped_SC_thresh(EigMeta200_all_grouped_SC_thresh<thresh) = 0;


meta_intdegseg_all = zeros(ngroups, 3, 50);

for i = 1:ngroups
    for j = 1:n_subjects

% % % % % With Threshold (80th Percentile)
        meta_intdegseg_all(i,1,j) = mean(EigMeta200_all_grouped_SC_thresh(i,int_range,j));
        meta_intdegseg_all(i,2,j) = mean(EigMeta200_all_grouped_SC_thresh(i,deg_range,j));
        meta_intdegseg_all(i,3,j) = mean(EigMeta200_all_grouped_SC_thresh(i,seg_range,j));

% % % % % % Without Threshold
%         meta_intdegseg_all(i,1,j) = mean(EigMeta200_all_grouped_SC(i,int_range,j));
%         meta_intdegseg_all(i,2,j) = mean(EigMeta200_all_grouped_SC(i,deg_range,j));
%         meta_intdegseg_all(i,3,j) = mean(EigMeta200_all_grouped_SC(i,seg_range,j));

    end
end

meta_intdegseg_all(isnan(meta_intdegseg_all)) = 0;
meta_intdegseg = mean(meta_intdegseg_all,3);

coords = meta_intdegseg./sum(meta_intdegseg,1);
coords = coords./sum(coords,2);
coords = coords(:,[1 3 2]);

% % % YMC Colormap
ycm_colors = zeros(size(coords));  % Initialize colors
ycm_colors(:, 1) = (1 - coords(:, 1)); % Cyan
ycm_colors(:, 2) = (1 - coords(:, 2)); % Yellow
ycm_colors(:, 3) = (1 - coords(:, 3)); % Magenta


coords = meta_intdegseg./sum(meta_intdegseg,2);

wlimits = ternary_axes_limits;

FIG = figure();
vgen  = { 'wlimits',       wlimits ,... % Axes will match wlimits ranges
'titlelabels', {'Int','Seg','Deg'}, ... % custom labels
'titlerotation', [0,0,0], ... % Set all titles to horizontal
'link_color', {'tick','title'},... % Link all axes colors
'titleshift',[ -0.18, 0, 0.18; 0.085, -0.11, 0.085 ]... % shift titles
'tickshift', [-0.02, -0.02, 0.0; -0.01,-0.00,0]
};
% Ternary Axes Outline   - Passed directly to plot3()
vout  = { 'LineWidth', 3, 'LineStyle', '-','Color','k'};
% Ternary gridlines  - Passed directly to plot3()
vgrid = { 'LineStyle','--','LineWidth',0.5, 'Color',[0 0 0 0.2] };
% Ternary Tick Line - Passed directly to plot3()
vtick_line = { 'LineStyle', '-', 'LineWidth', 1.5, 'Color',[0 0 0 1.0] };
% Ternary Tick Labels - Passed directly to text()
vtick_label = { 'FontWeight','Bold', 'FontSize', 8 };
% Ternary Axes Label  - Passed directly to text()
vlab  = { 'FontWeight','normal', 'FontSize', 14 };
h = ternary_axes( vgen, vout, vgrid, vtick_line, vtick_label, vlab );
[pH, ~, flat_coords] = ternary_scatter3( wlimits, 'l', coords(:,1), 'b', coords(:,3), coords(:,2));
pH.MarkerEdgeColor = [0, 0, 0];
pH.CData = ycm_colors;
pH.SizeData = 200;
colorbar off;

hold on;
for i = 1:length(domain_names_abbv)
    text(flat_coords(i,1)+0.015, flat_coords(i,2), flat_coords(i,3), domain_names_abbv(i))
end

set(FIG, 'Position',[1,49,1920,955]);

