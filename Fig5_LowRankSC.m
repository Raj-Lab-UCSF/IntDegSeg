%% Figure 5: Low Rank SC

scale = 200;

% disp(scale);
filename = sprintf('MICA_schaefer%d_SCFC_struct.mat',scale);
load(filename);

n_subjects = length(MICA);
nroi = scale+14;

SC_all = zeros(nroi, nroi, n_subjects);

GraphMetrics = struct('Full',[],'Int',[],'Deg',[],'Seg',[]);

for sub = 1:n_subjects
    disp(sub);
    SC = MICA(sub).SC;
    SC = SC./norm(SC,'fro');
    SC_all(:,:,sub) = SC;
    GraphMetrics(sub).Full = analyze_network(SC);
end
SC_consensus = mean(SC_all,3);
[~, SC_U_consensus] = graph_laplacian(SC_consensus,'normalized');

%%%% Determine from all subjects:
IntDegSeg_Consensus = zeros(3, nroi, n_subjects);
IntDegSeg_sub_assigment = zeros(3,nroi);

for sub = 1:n_subjects

    [~, U, ev] = graph_laplacian(SC_all(:,:,sub), 'normalized');

    [ev_changepts, ~, ~, ~] = find_ev_IntDegSeg(ev, '2deriv', 3, 10);

    IntDegSeg_sub_assigment = zeros(3,nroi);
    IntDegSeg_sub_assigment(1,1:(ev_changepts(1)-1)) = 1;
    IntDegSeg_sub_assigment(2,ev_changepts(1):(ev_changepts(2)-1)) = 1;
    IntDegSeg_sub_assigment(3,ev_changepts(2):nroi) = 1;

    [~, sigma] = match_eigs(U, SC_U_consensus);

    IntDegSeg_Consensus(:,:,sub) = IntDegSeg_sub_assigment(:, sigma);
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

filter_weights = ones((scale+14),1);
filter_weights(deg_range) = 2;
filter_weights(seg_range) = 3;

%
for sub = 1:n_subjects
    disp(sub);

    SC = SC_all(:,:,sub);
    [~, SC_U] = graph_laplacian(SC,'normalized');

    [~, sigma] = match_eigs(SC_U, SC_U_consensus);

    SC_IntDegSeg = filterHarmonics(SC, filter_weights, sigma);

    SC_Int = SC_IntDegSeg(:,:,1)./norm(SC_IntDegSeg(:,:,1),'fro');
    SC_Deg = SC_IntDegSeg(:,:,2)./norm(SC_IntDegSeg(:,:,2),'fro');
    SC_Seg = SC_IntDegSeg(:,:,3)./norm(SC_IntDegSeg(:,:,3),'fro');

    GraphMetrics(sub).Int = analyze_network(SC_Int);
    GraphMetrics(sub).Deg = analyze_network(SC_Deg);
    GraphMetrics(sub).Seg = analyze_network(SC_Seg);


end

% save(sprintf('PreComputed_Data/LowRankSC_Scale-%d_GraphMetrics_MICA.mat',scale),'GraphMetrics');

%%
% % Unpack Structure:
% load('PreComputed_Data/LowRankSC_Scale-200_GraphMetrics_MICA.mat');
n_subjects = length(GraphMetrics);
scale = 200;

SmallWorldness = zeros(n_subjects, 4);
RC_AUC = zeros(n_subjects, 4);
Eff = zeros(n_subjects, 4);
ModQ = zeros(n_subjects, 4);
Core = zeros(n_subjects, 4);


for sub = 1:n_subjects


    SmallWorldness(sub,1) = GraphMetrics(sub).Full.global.small_world.sigma;
    SmallWorldness(sub,2) = GraphMetrics(sub).Int.global.small_world.sigma;
    SmallWorldness(sub,3) = GraphMetrics(sub).Deg.global.small_world.sigma;
    SmallWorldness(sub,4) = GraphMetrics(sub).Seg.global.small_world.sigma;

    RC_Full = [GraphMetrics(sub).Full.global.rich_club_curve];
    RC_Int = [GraphMetrics(sub).Int.global.rich_club_curve];
    RC_Deg = [GraphMetrics(sub).Deg.global.rich_club_curve];
    RC_Seg = [GraphMetrics(sub).Seg.global.rich_club_curve];

    RC_Full(isnan(RC_Full)) = 0;
    RC_Int(isnan(RC_Int)) = 0;
    RC_Deg(isnan(RC_Deg)) = 0;
    RC_Seg(isnan(RC_Seg)) = 0;

    RC_AUC(sub,1) = trapz(RC_Full);
    RC_AUC(sub,2) = trapz(RC_Int);
    RC_AUC(sub,3) = trapz(RC_Deg);
    RC_AUC(sub,4) = trapz(RC_Seg);


    Eff(sub,1) = GraphMetrics(sub).Full.global.efficiency;
    Eff(sub,2) = GraphMetrics(sub).Int.global.efficiency;
    Eff(sub,3) = GraphMetrics(sub).Deg.global.efficiency;
    Eff(sub,4) = GraphMetrics(sub).Seg.global.efficiency;


    ModQ(sub,1) = GraphMetrics(sub).Full.global.modularity_Q;
    ModQ(sub,2) = GraphMetrics(sub).Int.global.modularity_Q;
    ModQ(sub,3) = GraphMetrics(sub).Deg.global.modularity_Q;
    ModQ(sub,4) = GraphMetrics(sub).Seg.global.modularity_Q;


    Core(sub,1) = GraphMetrics(sub).Full.global.coreness_q;
    Core(sub,2) = GraphMetrics(sub).Int.global.coreness_q;
    Core(sub,3) = GraphMetrics(sub).Deg.global.coreness_q;
    Core(sub,4) = GraphMetrics(sub).Seg.global.coreness_q;

end

names = {'Full','INT','DEG','SEG'};

full_color = [0.3, 0.3, 0.3];
int_color = [0.3, 0.9, 0.9];
deg_color = [0.9, 0.9, 0.3];
seg_color = [0.9, 0.3, 0.9];
colors = [full_color; int_color; deg_color; seg_color];




%% Violin: Modularity

measure_name = 'Modularity';

FIG = figure();
violinplot(ModQ, names, 'ViolinColor',colors);
title(measure_name);
axis square;
% set(gca, 'YScale', 'log');
set(FIG, 'Position',[1,49,1920,955]);

% % % Stats:
p = anova1(ModQ);

p_mat = zeros(4);
for i = 1:4
    for j = 1:4
        [~, p_mat(i,j)] = ttest2(ModQ(:,i), ModQ(:,j));
    end
end

%% Violin: Efficiency

measure_name = 'Efficiency';

FIG = figure();
violinplot(Eff, names, 'ViolinColor',colors);
title(measure_name);
axis square;
ylim([3.5E-3, 9E-3]);
set(FIG, 'Position',[1,49,1920,955]);

% % % Stats:
p = anova1(Eff);

p_mat = zeros(4);
for i = 1:4
    for j = 1:4
        [~, p_mat(i,j)] = ttest2(Eff(:,i), Eff(:,j));
    end
end

%% Violin: Small Worldness

measure_name = 'SmallWorldness';

FIG = figure();
violinplot(SmallWorldness, names, 'ViolinColor',colors);
title(measure_name);
axis square;
% set(gca, 'YScale', 'log');
set(FIG, 'Position',[1,49,1920,955]);

% % % Stats:
p = anova1(SmallWorldness);

p_mat = zeros(4);
for i = 1:4
    for j = 1:4
        [~, p_mat(i,j)] = ttest2(SmallWorldness(:,i), SmallWorldness(:,j));
    end
end

%% Violin: Coreness

measure_name = 'Coreness';

FIG = figure();
violinplot(Core, names, 'ViolinColor',colors);
title(measure_name);
axis square;
ylim([0.25,0.7])
set(FIG, 'Position',[1,49,1920,955]);

% % % Stats:
p = anova1(Core);

p_mat = zeros(4);
for i = 1:4
    for j = 1:4
        [~, p_mat(i,j)] = ttest2(Core(:,i), Core(:,j));
    end
end

%% Violin: Rich Club AUC

measure_name = 'RC_{AUC}';

FIG = figure();
violinplot(RC_AUC, names, 'ViolinColor',colors);
title(measure_name);
axis square;
set(gca, 'YScale', 'log');
set(FIG, 'Position',[1,49,1920,955]);

% % % Stats:
p = anova1(RC_AUC);

p_mat = zeros(4);
for i = 1:4
    for j = 1:4
        [~, p_mat(i,j)] = ttest2(RC_AUC(:,i), RC_AUC(:,j));
    end
end



%% Violin: Radius

measure_name = 'Radius';

FIG = figure();
violinplot(Radius, names, 'ViolinColor',colors);
title(measure_name);
axis square;
ylim([5E2, 1E6]);
set(gca, 'YScale', 'log');
set(FIG, 'Position',[1,49,1920,955]);



%% Violin: Assortativity

measure_name = 'Assortativity';

FIG = figure();
violinplot(Assort, names, 'ViolinColor',colors);
title(measure_name);
axis square;
% set(gca, 'YScale', 'log');
set(FIG, 'Position',[1,49,1920,955]);





%% Filtered SC

FIG = figure();
imagesc(SC_consensus);
clim([0, prctile(nonzeros(SC_consensus),98, 'all')]);
xticks([]);
yticks([]);
xticklabels([]);
yticklabels([]);
axis square;
set(FIG, 'Position',[1,49,1920,955])


FIG = figure();
imagesc(SC_IntDegSeg(:,:,1));
clim([0, prctile(nonzeros(SC_IntDegSeg(:,:,1)),98, 'all')]);
xticks([]);
yticks([]);
xticklabels([]);
yticklabels([]);
axis square;
set(FIG, 'Position',[1,49,1920,955])


FIG = figure();
imagesc(SC_IntDegSeg(:,:,2));
clim([0, prctile(nonzeros(SC_IntDegSeg(:,:,2)),98, 'all')]);
xticks([]);
yticks([]);
xticklabels([]);
yticklabels([]);
axis square;
set(FIG, 'Position',[1,49,1920,955])

FIG = figure();
imagesc(SC_IntDegSeg(:,:,3));
clim([0, prctile(nonzeros(SC_IntDegSeg(:,:,3)),98, 'all')]);
xticks([]);
yticks([]);
xticklabels([]);
yticklabels([]);
axis square;
set(FIG, 'Position',[1,49,1920,955])











