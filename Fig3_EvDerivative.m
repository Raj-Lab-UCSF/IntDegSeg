%% Figure Script for the Eigenvalue Derivative Analysis

%% Setup Data:

scale = 200;

filename = 'MICA_schaefer200_SCFC_struct.mat';

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
[SC_L, SC_U_consensus, SC_ev_consensus] = graph_laplacian(SC_consensus, 'normalized');


int_color = [0.1, 0.9, 0.9];
deg_color = [0.9, 0.9, 0.1];
seg_color = [0.9, 0.1, 0.9];

%% Find Int Deg Seg Regimes:

% %%%% Determine Int Deg Seg Regimes from Consensus Connectome eigenvalues:
% [ev_changepts, ev_fit, ev_deriv] = find_ev_IntDegSeg(SC_ev_consensus, '2deriv');
% 
% int_range = (1:(ev_changepts(1)-1))';
% deg_range = (ev_changepts(1):(ev_changepts(2)-1))';
% seg_range = (ev_changepts(2):nroi)';


%%%% Determine from all subjects:
IntDegSeg_Consensus = zeros(3, nroi, n_subjects);
IntDegSeg_sub_assigment = zeros(3,nroi);

for i = 1:n_subjects

    [~, U, ev] = graph_laplacian(SC_all(:,:,i), 'normalized');

    
    [ev_changepts, ~, ~, ~] = find_ev_IntDegSeg(ev, '2deriv', 3, 10);


    
    IntDegSeg_sub_assigment = zeros(3,nroi);
    IntDegSeg_sub_assigment(1,1:(ev_changepts(1)-1)) = 1;
    IntDegSeg_sub_assigment(2,ev_changepts(1):(ev_changepts(2)-1)) = 1;
    IntDegSeg_sub_assigment(3,ev_changepts(2):nroi) = 1;

    [~, sigma] = match_eigs(U, SC_U_consensus);

    IntDegSeg_Consensus(:,:,i) = IntDegSeg_sub_assigment(:, sigma);
end

IntDegSeg_Consensus_mean = mean(IntDegSeg_Consensus,3);

% % Find the consensus crossings:
pts=findchangepts(sign(IntDegSeg_Consensus_mean(2,:) - IntDegSeg_Consensus_mean(1,:)),'MaxNumChanges',3);
int_range_bound = min(pts);
pts=findchangepts(sign(IntDegSeg_Consensus_mean(3,:) - IntDegSeg_Consensus_mean(2,:)),'MaxNumChanges',3);
seg_range_bound = max(pts);

% % Define Regimes:
int_range = 1:(int_range_bound-1);
deg_range = int_range_bound:(seg_range_bound-1);
seg_range = seg_range_bound:nroi;


lw = 3; %linewidth

% % Show Regime Consensus plot with crossings:
FIG = figure();

t = tiledlayout(1,1);

ax1 = axes(t);
p_int = plot(ax1, 1:nroi, IntDegSeg_Consensus_mean(1,:)');
p_int.Color = int_color;
p_int.LineWidth = lw;
ax1.XLim = [0 scale+20];
ax1.Box = 'on';
set(ax1, 'XTickLabel',[]); set(ax1, 'YTickLabel',[]);

ax2 = axes(t);
p_deg = plot(ax2, 1:nroi, IntDegSeg_Consensus_mean(2,:)');
p_deg.Color = deg_color;
p_deg.LineWidth = lw;
ax2.XLim = [0 scale+20];
ax2.Color = 'none';
ax2.Box = 'off';
set(ax2, 'XTickLabel',[]); set(ax2, 'YTickLabel',[]);

ax3 = axes(t);
p_seg = plot(ax3, 1:nroi, IntDegSeg_Consensus_mean(3,:)');
p_seg.Color = seg_color;
p_seg.LineWidth = lw;
ax3.XLim = [0 scale+20];
ax3.Color = 'none';
ax3.Box = 'off';

xline([int_range_bound, seg_range_bound], '--k', 'LineWidth',3);

fig_width = 6;
fig_height = 2;

pbaspect(ax1, [fig_width, fig_height, 1]);
pbaspect(ax2, [fig_width, fig_height, 1]);
pbaspect(ax3, [fig_width, fig_height, 1]);



%% Pannel A: Eigenvalues and their Derivative:
fig_dir = 'C:\Users\flame\Box\Projects\Eigenmode_Variability\Figures\Fig4_EVDerivative\';

[~, ev_fit, ev_deriv, ev_2deriv] = find_ev_IntDegSeg(SC_ev_consensus, '2deriv', 3, 10);

sc_color = [1 0.25 0.5];
fit_color = [0.25 0.5 1];
deriv_color = [0.5 .2 .8];
freq_color = [0.1 0.4 0.1];
lw = 3;

FIG = figure();

ax1 = axes;
p1 = plot(ax1, 1:nroi, ev_deriv);
p1.Color = deriv_color;
p1.LineWidth = lw;
ax1.YAxisLocation = 'right';
ax1.XLim = [0 scale+20];
ax1.YColor = deriv_color;
ax1.Box = 'off';
ylabel('Gap-Spectrum');
set(ax1, 'XTickLabel',[]); 


ax2 = axes;
p2 = plot(ax2, 1:nroi, SC_ev_consensus);
p2.Color = sc_color;
p2.LineWidth = lw;
ax2.Color = 'none';
ax2.XLim = [0 scale+20];
ax2.YColor = sc_color;
ylabel('Eigenvalue');
xlabel('Eigenmode Index');
ax2.YLim = [0 1.5];
ax2.Box = 'off';
set(ax2, 'XTickLabel',[]);

xline([int_range_bound, seg_range_bound], '--k', 'LineWidth',3);



FIG = figure();
p3 = plot(1:nroi, ev_2deriv);
p3.Color = fit_color;
p3.LineWidth = lw;
ylabel('2nd Order Gap-Spectrum');
xlabel('Eigenmode Index');
title('Second-Order Gap-Spectrum');


%% % Find Zerocrossrate, entropy, min-cut-max-flow, and sparsity

load(sprintf('Schaefer%d_Adjacency.mat',scale));

Adja_norm = Adja./norm(Adja,'fro');
L_adja = graph_laplacian(Adja_norm,'normalized');


thresh = 1e-3;

SC_ev_roughness = zeros(nroi,n_subjects);
SC_ev_zeroX = zeros(nroi,n_subjects);
SC_ev_sparsity = zeros(nroi,n_subjects);
SC_ev_mincut = zeros(nroi,n_subjects);
ev_deriv_all = zeros(nroi,n_subjects);

for sub = 1:n_subjects

    SC = SC_all(:,:,sub);
    [~, SC_U, SC_ev] = graph_laplacian(SC, 'normalized');

    [SC_U_matched, sigma] = match_eigs(SC_U, SC_U_consensus);


    SC_ev_roughness(:,sub) = vecnorm(L_adja*SC_U_matched)';

    SC_ev_mincut(:,sub) = ev_zeroXings(SC, SC_U_matched, thresh);

    for i = 1:nroi
    
        evec = SC_U_matched(:,i);
        
        evec_thresh = evec;
        evec_thresh(abs(evec_thresh) <=thresh) = 0;
    
        SC_ev_zeroX(i,sub) = zerocrossrate(evec_thresh);
        
        evec_thresh(abs(evec_thresh) >thresh) = 1;
    
        SC_ev_sparsity(i,sub) = (nroi - sum(evec_thresh))/nroi;
    
    end

end


SC_ev_roughness_median = median(SC_ev_roughness,2);
SC_ev_mincut_median = median(SC_ev_mincut,2);
SC_ev_zeroX_median = median(SC_ev_zeroX,2);
SC_ev_sparsity_median = median(SC_ev_sparsity,2);

% % % %  Note: mean and median produce almost identical results
% SC_ev_entropy_median = mean(SC_ev_entropy,2);
% SC_ev_mincut_median = mean(SC_ev_mincut,2);
% SC_ev_zeroX_median = mean(SC_ev_zeroX,2);
% SC_ev_sparsity_median = mean(SC_ev_sparsity,2);



%% Linear Regression for all variables:


% % consensus Derivative:
X_int = [ones(size(int_range')), ev_deriv(int_range)]; % Design matrix
X_deg = [ones(size(deg_range')), ev_deriv(deg_range)]; % Design matrix
X_seg = [ones(size(seg_range')), ev_deriv(seg_range)]; % Design matrix



% % % %  Linear Regression Eigenvector Roughness:
b = X_int \ SC_ev_roughness_median(int_range); % Regression coefficients (b0, b1)
LRfit_roughness_int = X_int * b; % Fitted values

b = X_deg \ SC_ev_roughness_median(deg_range); % Regression coefficients (b0, b1)
LRfit_roughness_deg = X_deg * b; % Fitted values

b = X_seg \ SC_ev_roughness_median(seg_range); % Regression coefficients (b0, b1)
LRfit_roughness_seg = X_seg * b; % Fitted values


% % % %  Linear Regression Eigenvector Sparsity:
b = X_int \ SC_ev_sparsity_median(int_range); % Regression coefficients (b0, b1)
LRfit_sparsity_int = X_int * b; % Fitted values

b = X_deg \ SC_ev_sparsity_median(deg_range); % Regression coefficients (b0, b1)
LRfit_sparsity_deg = X_deg * b; % Fitted values

b = X_seg \ SC_ev_sparsity_median(seg_range); % Regression coefficients (b0, b1)
LRfit_sparsity_seg = X_seg * b; % Fitted values


% % % %  Linear Regression Eigenvector Min-Cut:
b = X_int \ SC_ev_mincut_median(int_range); % Regression coefficients (b0, b1)
LRfit_mincut_int = X_int * b; % Fitted values

b = X_deg \ SC_ev_mincut_median(deg_range); % Regression coefficients (b0, b1)
LRfit_mincut_deg = X_deg * b; % Fitted values

b = X_seg \ SC_ev_mincut_median(seg_range); % Regression coefficients (b0, b1)
LRfit_mincut_seg = X_seg * b; % Fitted values


% % % %  Linear Regression Eigenvector Zero Cross Rate:
b = X_int \ SC_ev_zeroX_median(int_range); % Regression coefficients (b0, b1)
LRfit_zeroX_int = X_int * b; % Fitted values

b = X_deg \ SC_ev_zeroX_median(deg_range); % Regression coefficients (b0, b1)
LRfit_zeroX_deg = X_deg * b; % Fitted values

b = X_seg \ SC_ev_zeroX_median(seg_range); % Regression coefficients (b0, b1)
LRfit_zeroX_seg = X_seg * b; % Fitted values


int_color = [0.3, 0.9, 0.9];
deg_color = [0.9, 0.9, 0.3];
seg_color = [0.9, 0.3, 0.9];
deriv_color = [0.5 0.2 0.8];
sz = 300;
lw = 3;

%% % % Panel: Eigenvector Roughness
FIG = figure();
hold on;
s1 = scatter(ev_deriv(int_range), SC_ev_roughness_median(int_range),  sz, 'filled', 'MarkerFaceColor', int_color, 'MarkerEdgeColor', 'black');
p1 = plot(ev_deriv(int_range), LRfit_roughness_int, ':','LineWidth',lw, 'Color',int_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Roughness');
ylabel('Roughness');
xlabel('Gap-Spectrum');

FIG = figure();
hold on;
s2 = scatter(ev_deriv(deg_range), SC_ev_roughness_median(deg_range),  sz, 'filled', 'MarkerFaceColor', deg_color, 'MarkerEdgeColor', 'black');
p2 = plot(ev_deriv(deg_range), LRfit_roughness_deg, ':','LineWidth',lw, 'Color',deg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Roughness');
ylabel('Roughness');
xlabel('Gap-Spectrum');

FIG = figure();
hold on;
s3 = scatter(ev_deriv(seg_range), SC_ev_roughness_median(seg_range),  sz, 'filled', 'MarkerFaceColor', seg_color, 'MarkerEdgeColor', 'black');
p3 = plot(ev_deriv(seg_range), LRfit_roughness_seg, ':','LineWidth',lw, 'Color',seg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Roughness');
ylabel('Roughness');
xlabel('Gap-Spectrum');


%% % % Panel: Eigenvector Sparsity
FIG = figure();
hold on;
s1 = scatter(ev_deriv(int_range), SC_ev_sparsity_median(int_range),  sz, 'filled', 'MarkerFaceColor', int_color, 'MarkerEdgeColor', 'black');
p1 = plot(ev_deriv(int_range), LRfit_sparsity_int, ':','LineWidth',lw, 'Color',int_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Sparsity');
ylabel('Sparsity');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);


FIG = figure();
hold on;
s2 = scatter(ev_deriv(deg_range), SC_ev_sparsity_median(deg_range),  sz, 'filled', 'MarkerFaceColor', deg_color, 'MarkerEdgeColor', 'black');
p2 = plot(ev_deriv(deg_range), LRfit_sparsity_deg, ':','LineWidth',lw, 'Color',deg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Sparsity');
ylabel('Sparsity');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);

FIG = figure();
hold on;
s3 = scatter(ev_deriv(seg_range), SC_ev_sparsity_median(seg_range),  sz, 'filled', 'MarkerFaceColor', seg_color, 'MarkerEdgeColor', 'black');
p3 = plot(ev_deriv(seg_range), LRfit_sparsity_seg, ':','LineWidth',lw, 'Color',seg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Sparsity');
ylabel('Sparsity');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);

%% % % Panel: Eigenvector MinCut

FIG = figure();
hold on;
s1 = scatter(ev_deriv(int_range), SC_ev_mincut_median(int_range),  sz, 'filled', 'MarkerFaceColor', int_color, 'MarkerEdgeColor', 'black');
p1 = plot(ev_deriv(int_range), LRfit_mincut_int, ':','LineWidth',lw, 'Color',int_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Min-Cut');
ylabel('Min-Cut');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);

FIG = figure();
hold on;
s2 = scatter(ev_deriv(deg_range), SC_ev_mincut_median(deg_range),  sz, 'filled', 'MarkerFaceColor', deg_color, 'MarkerEdgeColor', 'black');
p2 = plot(ev_deriv(deg_range), LRfit_mincut_deg, ':','LineWidth',lw, 'Color',deg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Min-Cut');
ylabel('Min-Cut');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);


FIG = figure();
hold on;
s3 = scatter(ev_deriv(seg_range), SC_ev_mincut_median(seg_range),  sz,  'filled', 'MarkerFaceColor', seg_color, 'MarkerEdgeColor', 'black');
p3 = plot(ev_deriv(seg_range), LRfit_mincut_seg, ':','LineWidth',lw, 'Color',seg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Min-Cut');
ylabel('Min-Cut');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);


%% % % Panel: Eigenvector Zero Cross Rate

FIG = figure();
hold on;
s1 = scatter(ev_deriv(int_range), SC_ev_zeroX_median(int_range),  sz, 'filled', 'MarkerFaceColor', int_color, 'MarkerEdgeColor', 'black');
p1 = plot(ev_deriv(int_range), LRfit_zeroX_int, ':','LineWidth',lw, 'Color',int_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Zero Cross Rate');
ylabel('Zero Cross Rate');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);

FIG = figure();
hold on;
s2 = scatter(ev_deriv(deg_range), SC_ev_zeroX_median(deg_range),  sz, 'filled', 'MarkerFaceColor', deg_color, 'MarkerEdgeColor', 'black');
p2 = plot(ev_deriv(deg_range), LRfit_zeroX_deg, ':','LineWidth',lw, 'Color',deg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Zero Cross Rate');
ylabel('Zero Cross Rate');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);

FIG = figure();
hold on;
s3 = scatter(ev_deriv(seg_range), SC_ev_zeroX_median(seg_range),  sz, 'filled', 'MarkerFaceColor', seg_color, 'MarkerEdgeColor', 'black');
p3 = plot(ev_deriv(seg_range), LRfit_zeroX_seg, ':','LineWidth',lw, 'Color',seg_color-0.3);
hold off;
set(gca,'XColor',deriv_color);
box on;
title('Eigenvector Zero Cross Rate');
ylabel('Zero Cross Rate');
xlabel('Gap-Spectrum');
pbaspect([fig_width, fig_height, 1]);
set(FIG, 'Position',[1,49,1920,955]);

%% Correlation Coefficients:

disp('Roughness:')
[r, p] = corr(ev_deriv(int_range), SC_ev_roughness_median(int_range))
[r, p] = corr(ev_deriv(deg_range), SC_ev_roughness_median(deg_range))
[r, p] = corr(ev_deriv(seg_range), SC_ev_roughness_median(seg_range))

disp('Sparsity:')
[r, p] = corr(ev_deriv(int_range), SC_ev_sparsity_median(int_range))
[r, p] = corr(ev_deriv(deg_range), SC_ev_sparsity_median(deg_range))
[r, p] = corr(ev_deriv(seg_range), SC_ev_sparsity_median(seg_range))

disp('MinCut:')
[r, p] = corr(ev_deriv(int_range), SC_ev_mincut_median(int_range))
[r, p] = corr(ev_deriv(deg_range), SC_ev_mincut_median(deg_range))
[r, p] = corr(ev_deriv(seg_range), SC_ev_mincut_median(seg_range))

disp('ZeroX:')
[r, p] = corr(ev_deriv(int_range), SC_ev_zeroX_median(int_range))
[r, p] = corr(ev_deriv(deg_range), SC_ev_zeroX_median(deg_range))
[r, p] = corr(ev_deriv(seg_range), SC_ev_zeroX_median(seg_range))

