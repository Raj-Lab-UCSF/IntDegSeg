%% Experiment 12: Noise Stability Ananlysis

%% Perturbation:

rng('shuffle','twister');


U_Scales_Noise = struct('scale',[],'noise_sigma',[],'UU_diag_SC_vs_SCnoisy_all',[],'all_combo_diagonals',[]);

d = 2.5e-4;
std_multiplier_list = 0:2.5e-4:0.05;
var_multiplier_list = std_multiplier_list.^2;
k_multipliers = length(std_multiplier_list);


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

SC_template = mean(SC_all,3);
[~, SC_U_template, SC_ev_template] = graph_laplacian(SC_template, 'normalized');

UU_diag_SC_vs_SCnoisy_all = ones(n_subjects, nroi, k_multipliers);



for k = 2:k_multipliers
    disp(k);
    r = makedist('Rician','s',0,'sigma',std_multiplier_list(k));
    noise_mat = random(r,[nroi, nroi]);
    noise_mat = (triu(noise_mat,1) + triu(noise_mat,1)').*~eye(nroi);

    for i = 1:n_subjects

        SC = SC_all(:,:,i);
        [L_SC, SC_U_sub, ~] = graph_laplacian(SC, 'normalized');
        SC_U_sub_matched = match_eigs(SC_U_sub, SC_U_template);

        SC_noisy = (SC + noise_mat);
        SC_noisy = SC_noisy./norm(SC_noisy, 'fro');
        [L_SC_noisy, SC_noisy_U_sub, ~] = graph_laplacian(SC_noisy, 'normalized');
        
        SC_noisy_U_sub_matched2Sub = match_eigs(SC_noisy_U_sub, SC_U_sub_matched);
    
        UU_diag_SC_vs_SCnoisy_all(i,:,k) = abs(diag(SC_U_sub_matched' * SC_noisy_U_sub_matched2Sub));

    end

    
end

% save('PreComputed_Data/Stability_Analysis_Results.mat','UU_diag_SC_vs_SCnoisy_all','std_multiplier_list','var_multiplier_list','k_multipliers' ,'SC_template','SC_U_template','SC_ev_template');


%% Similarity by regime

% load('PreComputed_Data/Stability_Analysis_Results.mat');

n_subjects = length(MICA);

% % From EV Regime Consensus:
int_range = 1:26;
deg_range = 27:187;
seg_range = 188:214;

int_color = [0.3, 0.9, 0.9];
deg_color = [0.9, 0.9, 0.3];
seg_color = [0.9, 0.3, 0.9];

mean_SC_weight = mean(nonzeros(triu(SC_template,1)));

k_multipliers = 10;

Perturb_Int = zeros(n_subjects, k_multipliers);
Perturb_Deg = zeros(n_subjects, k_multipliers);
Perturb_Seg = zeros(n_subjects, k_multipliers);

for k = 1:k_multipliers

    Perturb_Int(:,k) = mean(UU_diag_SC_vs_SCnoisy_all(:, int_range,k),2);
    Perturb_Deg(:,k) = mean(UU_diag_SC_vs_SCnoisy_all(:, deg_range,k),2);
    Perturb_Seg(:,k) = mean(UU_diag_SC_vs_SCnoisy_all(:, seg_range,k),2);

end


FIG = figure();
hold on;


plot_iqr(std_multiplier_list(1:k_multipliers), Perturb_Int, 'median', int_color, true, 0.8);
plot_iqr(std_multiplier_list(1:k_multipliers), Perturb_Deg, 'median', deg_color, true, 0.8);
plot_iqr(std_multiplier_list(1:k_multipliers), Perturb_Seg, 'median', seg_color, true, 0.8);


xline(mean_SC_weight, 'k--','LineWidth',2);

hold off;
box on;
axis square;

xlim([0, std_multiplier_list(k_multipliers)]);
title('Network Perturbation');
xlabel('Rician Sigma');
ylabel('Perturbation Similarity');
set(FIG, 'Position',[1,49,1920,955]);


%% Compute Stability

customFit = fittype('m*x + 1', 'independent', 'x', 'dependent', 'y');

options = fitoptions('Method', 'NonlinearLeastSquares', ...
                     'StartPoint', [0], ...  % initial guess for [a, b]
                     'Lower', [-Inf], ...    % lower bounds for [a, b]
                     'Upper', [0]);      % upper bounds for [a, b]

k = 9; %67 = 1/3*201
x = std_multiplier_list(1:k)';

max_dPS_dP_lin = zeros(n_subjects,214);

for sub = 1:n_subjects
    for i = 1:nroi

        F0 = fit(x, log(squeeze(UU_diag_SC_vs_SCnoisy_all(sub,i,1:k))), customFit, options);
        max_dPS_dP_lin(sub,i) = F0.m;

    end
end


%% Visualize Stability Across Harmonics:

[~, ~, ev_deriv] = find_ev_IntDegSeg(SC_ev_template,'1deriv');

% % From EV Regime Consensus:
int_range = 1:26;
deg_range = 27:187;
seg_range = 188:214;

int_color = [0.3, 0.9, 0.9];
deg_color = [0.9, 0.9, 0.3];
seg_color = [0.9, 0.3, 0.9];
deriv_color = [0.5 .2 .8];

lw = 2;
X = 1:214;


FIG = figure(); 
t = tiledlayout(1,1);

ax1 = axes(t);
hold on;
p1 = plot_iqr(X(int_range),max_dPS_dP_lin(:,int_range), 'median', int_color, true,0.8);
p2 = plot_iqr(X(deg_range),max_dPS_dP_lin(:,deg_range), 'median', deg_color, true,0.8);
p3 = plot_iqr(X(seg_range),max_dPS_dP_lin(:,seg_range), 'median', seg_color, true,0.8);

p1.LineWidth = lw;
p2.LineWidth = lw;
p3.LineWidth = lw;
xline(int_range(end)+0.5, 'k--', 'LineWidth',lw);
xline(deg_range(end)+0.5, 'k--', 'LineWidth',lw);
ax1.Box = 'off';
ax1.XLim = [0 220];
ax1.YLim = [-1800 -600];
axis square;
hold off;


ax2 = axes(t);
p4 = plot(ax2, 1:nroi, ev_deriv);
p4.Color = deriv_color;
p4.LineWidth = 3;
ax2.YAxisLocation = 'right';
ax2.YColor = deriv_color;
ax2.Box = 'off';
ax2.Color = 'none';
set(ax2, 'XTickLabel',[]); 
ax2.XLim = [0 220];
axis square;
set(FIG, 'Position',[1,49,1920,955]);
ax1.YLabel = ylabel('Perturbation Derivative Maximum');
ax2.YLabel = ylabel('Eigenvalue Derivative');
title('Stability of Harmonics');
xlabel({'','Eigenmode Index'});


%% Stability vs Gap-Spectrum Histogram and Scatter Plot

r = zeros(n_subjects,1); p = r; sz = 20;
r_int = zeros(n_subjects,1); p_int = r; 
r_deg = zeros(n_subjects,1); p_deg = r; 
r_seg = zeros(n_subjects,1); p_seg = r; 

FIG = figure(); hold on;
for i = 1:n_subjects
    [r(i), p(i)] = corr(ev_deriv, max_dPS_dP_lin(i,:)');

    [r_int(i), p_int(i)] = corr(ev_deriv(int_range), max_dPS_dP_lin(i,int_range)');
    [r_deg(i), p_deg(i)] = corr(ev_deriv(deg_range), max_dPS_dP_lin(i,deg_range)');
    [r_seg(i), p_seg(i)] = corr(ev_deriv(seg_range), max_dPS_dP_lin(i,seg_range)');

    scatter(ev_deriv(int_range), max_dPS_dP_lin(i,int_range)', 30,'filled', 'MarkerFaceColor', int_color, 'MarkerEdgeColor','k');
    scatter(ev_deriv(deg_range), max_dPS_dP_lin(i,deg_range)', 30,'filled', 'MarkerFaceColor', deg_color, 'MarkerEdgeColor','k');
    scatter(ev_deriv(seg_range), max_dPS_dP_lin(i,seg_range)', 30,'filled', 'MarkerFaceColor', seg_color, 'MarkerEdgeColor','k');

end
% % Template Derivative:
ev_deriv_rep = repmat(ev_deriv,n_subjects,1);
X = [ones(size(ev_deriv_rep)), ev_deriv_rep]; % Design matrix

max_dPS_dP_transpose = max_dPS_dP_lin';
b = X \ max_dPS_dP_transpose(:); % Regression coefficients (b0, b1)
LRfit = X * b; % Fitted values
p1 = plot(ev_deriv_rep, LRfit, 'k--','LineWidth',lw);
hold off;
axis square;
box on;
title('All Subjects dLambda vs dPS');
xlabel('EV Deriv');
ylabel('dPS/dP');
set(FIG, 'Position',[1,49,1920,955]);
ylim([-1800, -600]);


%
% % % Scatter of just the median:
sz = 200;
median_max_dPS_dP = median(max_dPS_dP_lin)';
% median_max_dPS_dP = median(max_dPS_dP)';
[r_med, p_med] = corr(ev_deriv,  median_max_dPS_dP);

cmap = summer(3);

FIG = figure();
h = histogram(r);
h = histogram(r, 0.575:0.025:0.75);
ylim([0,25]);
p5 = xline(mean(r), 'k--');%,sprintf('mean=%0.2f', mean(r)));
p5.LineWidth = lw;
h.FaceColor = cmap(2,:);
axis square;
set(FIG, 'Position',[1,49,1920,955]);
title('Histogram of individual correlations');
xlabel("Pearson's Correlation");
ylabel('Subject Count');



%% Some extra plots...

% % Plots for median stability per regime:
% 
% FIG = figure(); hold on;
% X = [ones(length(int_range),1), ev_deriv(int_range)];
% b = X \ median_max_dPS_dP(int_range);
% LRfit = X * b; % Fitted values
% p1 = plot(ev_deriv(int_range), LRfit, ':','Color', [0 0 0],'LineWidth',lw);
% scatter(ev_deriv(int_range),median_max_dPS_dP(int_range), sz,'filled', 'MarkerFaceColor',int_color,'MarkerEdgeColor','k');
% hold off;
% axis square;
% box on;
% title('INT Median dLambda vs dPS');
% xlabel('EV Deriv');
% ylabel('dPS/dP');
% set(FIG, 'Position',[1,49,1920,955]);
% % ylim([-500, 0]);
% 
% FIG = figure(); hold on;
% X = [ones(length(deg_range),1), ev_deriv(deg_range)];
% b = X \ median_max_dPS_dP(deg_range);
% LRfit = X * b; % Fitted values
% p1 = plot(ev_deriv(deg_range), LRfit, ':','Color', [0 0 0],'LineWidth',lw);
% scatter(ev_deriv(deg_range),median_max_dPS_dP(deg_range), sz,'filled', 'MarkerFaceColor',deg_color,'MarkerEdgeColor','k');
% hold off;
% axis square;
% box on;
% title('DEG Median dLambda vs dPS');
% xlabel('EV Deriv');
% ylabel('dPS/dP');
% set(FIG, 'Position',[1,49,1920,955]);
% % ylim([-3000, 0]);
% 
% FIG = figure(); hold on;
% X = [ones(length(seg_range),1), ev_deriv(seg_range)];
% b = X \ median_max_dPS_dP(seg_range);
% LRfit = X * b; % Fitted values
% p1 = plot(ev_deriv(seg_range), LRfit, ':','Color', [0 0 0],'LineWidth',lw);
% scatter(ev_deriv(seg_range), median_max_dPS_dP(seg_range), sz,'filled', 'MarkerFaceColor',seg_color,'MarkerEdgeColor','k');
% hold off;
% axis square;
% box on;
% title('SEG Median dLambda vs dPS');
% xlabel('EV Deriv');
% ylabel('dPS/dP');
% set(FIG, 'Position',[1,49,1920,955]);
% % ylim([-2000, 0]);



% % Violin Plots per subject-wise harmonic stability vs regime
% data2plot = [r_int, r_deg, r_seg];
% 
% FIG = figure();
% violinplot(data2plot,{'INT','DEG','SEG'}, 'ViolinColor',[int_color; deg_color; seg_color]);
% 
% box on;
% 
% title('Subject Correlations Perturbation vs Regime');
% ylabel("Pearson's correlation");
% ylim([0 1]);
% fig_width = 3.75; fig_height = 1.45;
% pbaspect([fig_width, fig_height, 1]);
% set(FIG, 'Position',[1,49,1920,955]);

