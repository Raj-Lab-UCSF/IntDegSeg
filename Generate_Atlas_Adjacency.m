%% Atlas Adjacency Matrix

for s = 1:10

scale = s*100;
disp(scale);
nroi = scale+14;

atlas_path = 'C:\Users\flame\Box\atlases\Schaefer100_to_1000\';
% atlas_name = sprintf('Schaefer2018_%dParcels_7Networks_order_FSLMNI152_2mm.nii.gz', scale);
atlas_name = sprintf('Schaefer2018_%dParcels_7Networks_order_FSLMNI152_1mm.nii.gz', scale);

% 
nii = load_untouch_nii([atlas_path, atlas_name]);
atlas = nii.img;

img_sz = size(atlas);

Adja_cortex = zeros(scale,scale);

% % % Left Hemisphere
for i = 1:((scale/2)-1)

    ROI_A = zeros(img_sz);

    ROI_A(atlas == i) = 1;

    ROI_A_dist = bwdist(ROI_A, 'euclidean');

    for j = (i+1):(scale/2)

        ROI_B_idx = find(atlas == j);

        ROI_B_dist2A = ROI_A_dist(ROI_B_idx);

        ROI_B_dist2A(ROI_B_dist2A>=2) = 0;

        ROI_B_adja = (1./ROI_B_dist2A).^3;
        ROI_B_adja(isinf(ROI_B_adja)) = 0;
        ROI_B_adja(isnan(ROI_B_adja)) = 0;

        Adja_cortex(i,j) = sum(ROI_B_adja);
        Adja_cortex(j,i) = Adja_cortex(i,j);

    end
end

% % % Right Hemisphere
for i = ((scale/2)+1):(scale-1)

    ROI_A = zeros(img_sz);

    ROI_A(atlas == i) = 1;

    ROI_A_dist = bwdist(ROI_A, 'euclidean');

    for j = (i+1):scale

        ROI_B_idx = find(atlas == j);

        ROI_B_dist2A = ROI_A_dist(ROI_B_idx);

        ROI_B_dist2A(ROI_B_dist2A>=2) = 0;

        ROI_B_adja = (1./ROI_B_dist2A).^3;
        ROI_B_adja(isinf(ROI_B_adja)) = 0;
        ROI_B_adja(isnan(ROI_B_adja)) = 0;

        Adja_cortex(i,j) = sum(ROI_B_adja);
        Adja_cortex(j,i) = Adja_cortex(i,j);

    end
end

Adja = zeros(nroi,nroi);
Adja(15:end, 15:end) = Adja_cortex;

save(sprintf('Schaefer%d_Adjacency.mat', scale), 'Adja');

end


%% Generate Euclidean Distance from coordinates:

scale = 1000;
coords = cat(1, subcortex_coords, schaefer1000_coords);
savename = sprintf('Schaefer%d_EuclideanDistance.mat',scale);

combos = nchoosek(1:(scale+14), 2);

EucDist = zeros(scale+14);

for i = 1:length(combos)

    a = combos(i,1);
    b = combos(i,2);

    coord_a = coords(a,:)';
    coord_b = coords(b,:)';

    EucDist(a,b) = norm(coord_a - coord_b);
    EucDist(b,a) = EucDist(a,b);


end

save(savename, 'EucDist');



