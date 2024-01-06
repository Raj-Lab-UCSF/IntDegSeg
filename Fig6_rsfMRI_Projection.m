%% Script to visualize functional projection of eigenmode regimes

% Necessary for Ternary Plot:
addpath ternary_plots/data_plots/
addpath ternary_plots/axes_creation/
addpath ternary_plots/utilities/

%% Projection Analysis:

scale = 200;

filename = sprintf('MICA_schaefer%d_SCFC_struct.mat',scale);
load(filename);

n_subjects = length(MICA);
nroi = length(MICA(1).SC);
nt = length(MICA(end).fMRI_TS);

SC_all = zeros(nroi,nroi, n_subjects);

% % From EV Regime Consensus:
int_range = 1:26;
deg_range = 27:187;
seg_range = 188:214;

int_color = [0.3, 0.9, 0.9];
deg_color = [0.9, 0.9, 0.3];
seg_color = [0.9, 0.3, 0.9];

for i = 1:n_subjects
    SC = MICA(i).SC;
    SC_all(:,:,i) = SC./norm(SC,'fro');
end

SC_template = mean(SC_all,3);
[~, SC_U_template, SC_ev_template] = graph_laplacian(SC_template, 'normalized');

Proj_BORN = zeros(nt, nroi, n_subjects);
Proj_BORN_Energy = zeros(nroi, n_subjects);

Regional_Energy_BORN_INT = zeros(nroi, n_subjects);
Regional_Energy_BORN_DEG = zeros(nroi, n_subjects);
Regional_Energy_BORN_SEG = zeros(nroi, n_subjects);

TS_INT_BORN = zeros(nt, n_subjects);
TS_DEG_BORN = zeros(nt, n_subjects);
TS_SEG_BORN = zeros(nt, n_subjects);


for i = 1:n_subjects
    
    TS = MICA(i).fMRI_TS(1:nt,:); %note: TS is time x region

    SC = SC_all(:,:,i);
    [~, SC_U] = graph_laplacian(SC, 'normalized');
    SC_U_matched = match_eigs(SC_U, SC_U_template);

    SC_U_matched_born = abs(SC_U_matched).^2;
    SC_U_matched_born_norm = SC_U_matched_born./vecnorm(SC_U_matched_born);

    SC_U_INT = SC_U_matched_born_norm(:, int_range);
    SC_U_DEG = SC_U_matched_born_norm(:, deg_range);
    SC_U_SEG = SC_U_matched_born_norm(:, seg_range);

% % % % Creates Timeseries of the GFW projections averaged across a domain
    
    Proj_INT = TS*SC_U_INT;
    TS_INT_BORN(:,i) = mean(Proj_INT,2);
    Proj_DEG = TS*SC_U_DEG;
    TS_DEG_BORN(:,i) = mean(Proj_DEG,2);
    Proj_SEG = TS*SC_U_SEG;
    TS_SEG_BORN(:,i) = mean(Proj_SEG,2);


    Regional_Energy_BORN_INT(:,i) = (vecnorm(Proj_INT)*SC_U_INT')';
    Regional_Energy_BORN_DEG(:,i) = (vecnorm(Proj_DEG)*SC_U_DEG')';
    Regional_Energy_BORN_SEG(:,i) = (vecnorm(Proj_SEG)*SC_U_SEG')';


    Proj_BORN(:,:,i) = TS*SC_U_matched_born_norm;
    Proj_BORN_Energy(:,i) = vecnorm(Proj_BORN(:,:,i))';

end

Regional_Energy_BORN_INT_mean = mean(Regional_Energy_BORN_INT./length(int_range), 2);
Regional_Energy_BORN_DEG_mean = mean(Regional_Energy_BORN_DEG./length(deg_range), 2);
Regional_Energy_BORN_SEG_mean = mean(Regional_Energy_BORN_SEG./length(seg_range), 2);

%% GFW Timeseries

sub = 50;
lw = 1.5;
fig_width = 12;
fig_height = 1;

FIG = figure();
hold on;
p1 = plot(TS_INT_BORN(:, sub));
p1.Color = int_color;
p1.LineWidth = lw;

p2 = plot(TS_DEG_BORN(:, sub));
p2.Color = deg_color;
p2.LineWidth = lw;

p3 = plot(TS_SEG_BORN(:, sub));
p3.Color = seg_color;
p3.LineWidth = lw;

hold off;

box on;

title('Example Subject GFW Born Projection');
ylabel('Projection Magnitude');
xlabel('Timepoint');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);


%% Mean GFW Energy Histogram:

lw = 6;
fig_width = 3.61;
fig_height = 2.71;


% Mean Regime Energy
Proj_BORN_Energy_INT = mean(Proj_BORN_Energy(int_range,:));
Proj_BORN_Energy_DEG = mean(Proj_BORN_Energy(deg_range,:));
Proj_BORN_Energy_SEG = mean(Proj_BORN_Energy(seg_range,:));

data2plot = [Proj_BORN_Energy_INT', Proj_BORN_Energy_DEG', Proj_BORN_Energy_SEG'];

FIG = figure();
violinplot(data2plot,{'INT','DEG','SEG'}, 'ViolinColor',[int_color; deg_color; seg_color]);

box on;

title('Born GFW Energy: All Subjects, Mean Harmonic Energy');
ylabel('Projection Energy');
% ylabel('Population Density');

pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);


%% GFW Energy Across All Harmonics

fig_width = 4.65;
fig_height =0.85;

FIG = figure();
plot_iqr(int_range, Proj_BORN_Energy(int_range,:)', 'median',int_color, true,0.8);
plot_iqr(deg_range, Proj_BORN_Energy(deg_range,:)', 'median',deg_color, true,0.8);
plot_iqr(seg_range, Proj_BORN_Energy(seg_range,:)', 'median',seg_color, true,0.8);
title('Born fMRI Projection Energy');
xlabel('Harmonic Index');
ylabel('Harmonic Energy');
xline([int_range(end)+0.5, deg_range(end)+0.5], 'k--', 'LineWidth',3);
xlim([0, 214]);
pbaspect([fig_width, fig_height, 1]);
box on;
set(FIG, 'Position',[308,425,1203,279]);



%% Ternary Point Cloud
TS_regimes = [];
TS_regime_phasemap = [];
for sub = 1:50
% % % Born: Regime Mean then Abs
    TS_regimes_sub = [abs(mean(Proj_BORN(:,int_range,sub),2)), abs(mean(Proj_BORN(:,deg_range,sub),2)), abs(mean(Proj_BORN(:,seg_range,sub),2))];

    TS_regimes_sub_rel = TS_regimes_sub./sum(TS_regimes_sub,2);
    TS_regime_phasemap = cat(1, TS_regime_phasemap, TS_regimes_sub_rel);
end

coords = TS_regime_phasemap;

ycm_colors = zeros(size(coords));  % Initialize colors
ycm_colors(:, 1) = (1 - coords(:, 1)); % Cyan
ycm_colors(:, 2) = (1 - coords(:, 2)); % Yellow
ycm_colors(:, 3) = (1 - coords(:, 3)); % Magenta


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
pH = ternary_scatter3( wlimits, 'l', coords(:,1), 'b', coords(:,3), coords(:,2));
pH.MarkerEdgeColor = [0, 0, 0];
ycm_colors_ternary = ycm_colors(:, [1, 3, 2]);
pH.CData = ycm_colors_ternary;
pH.SizeData = 10;
colorbar off;

set(FIG, 'Position',[1,49,1920,955]);


%% Ternary PCA:

c = pca(coords);

FIG=figure(); h=bar(c(:,1)); axis square; h.FaceColor = 'flat'; h.CData=[int_color; deg_color; seg_color];
title('Phase Plot PC1 Loadings'); xticklabels({'INT','DEG','SEG'});

%% Marginal Plots:

bins = linspace(0,1,50);

FIG = figure();
histogram(coords(:,1), bins, 'FaceColor',int_color);
xlim([0,1]);
ylim([0,1e4]);
set(FIG, 'Position',[0,483,1920,221]);

FIG = figure();
histogram(coords(:,2), bins, 'FaceColor',deg_color);
xlim([0,1]);
ylim([0,1e4]);
set(FIG, 'Position',[0,483,1920,221]);

FIG = figure();
histogram(coords(:,3), bins, 'FaceColor',seg_color);
xlim([0,1]);
ylim([0,1e4]);
set(FIG, 'Position',[0,483,1920,221]);


%% Entropy Analysis:

addpath SampEn/

H_INT = zeros(nt, n_subjects);
H_DEG = zeros(nt, n_subjects);
H_SEG = zeros(nt, n_subjects);

SampEn_INT = zeros(n_subjects,1);
SampEn_DEG = zeros(n_subjects,1);
SampEn_SEG = zeros(n_subjects,1);

SpectralEn_INT = zeros(n_subjects,1);
SpectralEn_DEG = zeros(n_subjects,1);
SpectralEn_SEG = zeros(n_subjects,1);

HFD_INT = zeros(n_subjects,1);
HFD_DEG = zeros(n_subjects,1);
HFD_SEG = zeros(n_subjects,1);

% Params for Sample Entropy:
m = 5; r = 0.2;

for sub = 1:50

% % % BORN: Ratio to Maximum Entropy
    H_INT(:,sub) = presence_entropy(abs(Proj_BORN(:,int_range,sub)),'bits')./log2(length(int_range));
%     H_INT(:,sub) = presence_entropy(abs(Proj_BORN(:,int_range(2:end),sub)))./log(length(int_range(2:end)));
    H_DEG(:,sub) = presence_entropy(abs(Proj_BORN(:,deg_range,sub)),'bits')./log2(length(deg_range));
    H_SEG(:,sub) = presence_entropy(abs(Proj_BORN(:,seg_range,sub)),'bits')./log2(length(seg_range));



% % % BORN: Sample Entropy
    SampEn_tmp = zeros(nroi,1);
    for j = 1:nroi
        GFW_ts = Proj_BORN(:,j,sub)';
        SampEn_tmp(j) = sampen(GFW_ts, m, r,'chebychev');
% % % % Note: sampen has units in nats, 
% % % % in our article, we changed log > log2 to obtain units in bits.
    end

    SampEn_INT(sub) = mean(SampEn_tmp(int_range), 'omitnan');
    SampEn_DEG(sub) = mean(SampEn_tmp(deg_range), 'omitnan');
    SampEn_SEG(sub) = mean(SampEn_tmp(seg_range), 'omitnan');

end

lw = 6;
fig_width = 3.61;
fig_height = 2.71;

% % Ratio to Maximum Entropy Figure
data2plot = [mean(H_INT,1)', mean(H_DEG,1)', mean(H_SEG,1)'];

FIG = figure();
violinplot(data2plot,{'INT','DEG','SEG'}, 'ViolinColor',[int_color; deg_color; seg_color]);
box on;
title('Harmonic Spatial Diversity');
ylabel('Spatial Entropy (Ratio to Maximum Entropy)');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);


% % Sample Entropy Figure:
data2plot = [SampEn_INT, SampEn_DEG, SampEn_SEG];

FIG = figure();
violinplot(data2plot,{'INT','DEG','SEG'}, 'ViolinColor',[int_color; deg_color; seg_color]);
box on;
title(['Harmonic Temporal Diversity (m = ' num2str(m) '; r = ' num2str(r) ')']);
ylabel('Sample Entropy');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);

