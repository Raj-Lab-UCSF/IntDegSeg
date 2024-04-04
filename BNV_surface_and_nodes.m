function BNV_surface_and_nodes(data, NIFTI, node_coords, surface_node_idx, Config, surf, savename)
%BNV_surface_and_nodes Create a BNV image that colors both a brain surface
%and nodes (likely subcortex nodes) simultaneously.
%   Usage:
%   BNV_surface_and_nodes(data, NIFTI, node_coords, surface_node_idx, Config, surf, savename)
% 
%   Inputs:
%       data: a Nx1 vector of data where N is the number of nodes.
%       NIFTI: a nifti image file containing an atlas with 1:K regions to color.
%       node_coords: an Mx3 matrix with coordinates for the M regions to
%       surface_node_idx: a Nx1 vector where the surface regions are 1, and
%           the bubble node regions are 0.
%       config_file: Contains information about the View and Colormaping to
%       be used in BNV
%       surf: Options: 'ICMB152' or 'Ch2'
%       savename: Base name to save output and intermediate files
%       
%       Note: N == K + M
% 
%       Dependencies: NITFI toolbox, BrainNetViewer, ColorBrewer

nii = load_untouch_nii(NIFTI);
img = nii.img;



max_img = max(img,[],'all');
q1 = unique(img);
q2 = 0:max_img;

if ~isequal(q1,q2')
    error('This atlas is not colored in an ordinal way and/or is missing labels. Reconfigure the atlas so that it contains values equal to 0:K')
end

surface_idx = find(surface_node_idx);

if ~isequal(length(surface_idx), length(surface_node_idx))

    node_idx = find(~surface_node_idx);
    
    if ~isequal(length(node_idx), size(node_coords,1))
        error('The number of specified nodes (for BNV bubbles) is not equal to the number of input node coordinates.')
    end

    colored_atlas = zeros(size(img));

    for i = 1:length(surface_idx)
        colored_atlas(img==i) = data(surface_idx(i));
    end
    
    nii.img = colored_atlas;
    
    namesplit = strsplit(NIFTI,'.');
    colored_atlas_file = [namesplit{1}, '_', savename, '.nii.gz'];
    
    save_untouch_nii(nii, colored_atlas_file);
    
    
    nodes = [node_coords, data(node_idx), ones(length(node_idx),1), nan(length(node_idx),1)];
    
    node_file = [namesplit{1}, '_', savename, '.node' ];
    dlmwrite(node_file, nodes, 'delimiter','\t');
    
    outfile = [savename, '.tif'];
    
    load(Config);

    EC.edg.draw = 0;
    
    if isreal(sqrt(data)) % i.e. if there are only positive values in data...
        EC.vol.pn = min(nonzeros(data));
        EC.vol.px = max(data)+eps;
        EC.nod.color_map_low = 0;
        EC.nod.color_map_high = max(data)+eps;
        EC.nod.CM = EC.vol.CM;
        EC.nod.CMt = EC.vol.CM;
    else
        EC.vol.pn = min(abs(nonzeros(data)));
        EC.vol.px = max(data)+eps;
        EC.vol.nn = -1*min(abs(nonzeros(data)));
        EC.vol.nx = min(data)-eps;
        EC.nod.color_map_low = min(data)-eps;
        EC.nod.color_map_high = max(data)+eps;
        EC.vol.adjustCM = 0;

        if sum(data(node_idx)) == 0
            EC.nod.CM = [linspace(0,1,64)', linspace(0,1,64)',linspace(0,1,64)'];
        else
            g = [0.5, 0.5, 0.5];
            c = EC.nod.CM;
            c(30:34,:) = [g; g; g; g; g];
            EC.nod.CM = c;
        end
    end
    
    new_config = [savename, '_', Config];
    save(new_config, 'EC');
    BrainNet_MapCfg(['BrainMesh_' surf '_smoothed.nv'], colored_atlas_file, node_file, new_config, outfile);

else

    if ~isequal(length(surface_idx), max_img)
        error('The number of specified surface indices is not equal to the number of input data rois coordinates.')
    end

    colored_atlas = zeros(size(img));
    
    for i = 1:length(surface_idx)
        colored_atlas(img==i) = data(surface_idx(i));
    end
    
    nii.img = colored_atlas;
    
    namesplit = strsplit(NIFTI,'.');
    colored_atlas_file = [namesplit{1}, '_', savename, '.nii.gz'];
    
    save_untouch_nii(nii, colored_atlas_file);
    

    outfile = [savename, '.tif'];
    
    load(Config);

    EC.edg.draw = 0;
    
    if isreal(sqrt(data)) % i.e. if there are only positive values in data...
        EC.vol.pn = min(data);
        EC.vol.px = max(data)+eps;
    else
        EC.vol.pn = 0+eps;
        EC.vol.px = max(data)+eps;
        EC.vol.nn = 0-eps;
        EC.vol.nx = min(data)-eps;
        EC.vol.adjustCM = 0;
    end
    
    new_config = [savename, '_', Config];
    save(new_config, 'EC');
    BrainNet_MapCfg(['BrainMesh_' surf '_smoothed.nv'], colored_atlas_file, new_config, outfile);


end



end