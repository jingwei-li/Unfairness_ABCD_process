function ABCD_chord_learnedBBA(learned_BBA, outname)

% ABCD_chord_learnedBBA(learned_BBA, outname)
%
% 

load(learned_BBA, 'avg_learned_cov')

%% load Schaefer parcellation, group ROI-to-ROI BBA matrix into network-to-network BBA matrix
lh_annot = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage6', 'label', ...
    'lh.Schaefer2018_400Parcels_17Networks_order.annot');
rh_annot = fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', ...
    'Schaefer2018_LocalGlobal', 'Parcellations', 'FreeSurfer5.3', 'fsaverage6', 'label', ...
    'rh.Schaefer2018_400Parcels_17Networks_order.annot');
Schaefer_17 = {'VisCent', 'VisPeri', 'SomMotA', 'SomMotB', 'DorsAttnA', 'DorsAttnB', ...
    'SalVentAttnA', 'SalVentAttnB', 'LimbicA', 'LimbicB', 'ContC', 'ContA', 'ContB', ...
    'TempPar', 'DefaultC', 'DefaultA', 'DefaultB'};

[~, lh_c] = CBIG_read_annotation(lh_annot);
[~, rh_c] = CBIG_read_annotation(rh_annot);
ROI_per_net = cell(length(Schaefer_17), 1);
for n = 1:length(Schaefer_17)
    lh_idx = cellfun(@any, strfind(lh_c.struct_names, Schaefer_17{n}));
    rh_idx = cellfun(@any, strfind(rh_c.struct_names, Schaefer_17{n}));
    ROI_per_net{n} = [find(lh_idx)-1; find(rh_idx)-1+200];
end

learnedBBA = zeros(length(Schaefer_17)+1);
for n1 = 1:length(Schaefer_17)
    for n2 = 1:length(Schaefer_17)
        learnedBBA(n1, n2) = mean(mean(avg_learned_cov(ROI_per_net{n1}, ROI_per_net{n2})));
        learnedBBA(length(Schaefer_17)+1, n2) = mean(mean(avg_learned_cov(401:419, ROI_per_net{n2})));
    end
    learnedBBA(n1, length(Schaefer_17)+1) = mean(mean(avg_learned_cov(ROI_per_net{n1}, 401:419)));
end

%% load color table and network names
load(fullfile(getenv('CBIG_CODE_DIR'), 'stable_projects', 'brain_parcellation', ...
    'Yeo2011_fcMRI_clustering', '1000subjects_reference', '17NetworksColors.mat'))
colors(colors==0) = 1e-15;
networks17 = {'Visual A', 'Visual B', 'SomMot A', 'SomMot B', 'DorsAttn A', 'DorsAttn B', ...
    'Sal/VentAttn A', 'Sal/VentAttn B', 'Limbic A', 'Limbic B', 'Control C', 'Control A', ...
    'Control B', 'TempPar', 'Default C', 'Default A', 'Default B'};

f = figure;
x = learnedBBA; x(x<0) = 0;
circularGraph(x, 'Colormap', colors([2:end 1], :)./255, 'Label', [networks17 'Subcortex'])
pos = get(gcf, 'Position');
%set(gcf, 'Position', [0 0 pos(3:4)*2])
export_fig([outname '_pos'], '-png', '-nofontswap', '-a1');
set(gcf, 'color', 'w')
hgexport(f, [outname '_pos'])
close
    
f = figure;
x = -learnedBBA; x(x<0) = 0;
circularGraph(x, 'Colormap', colors([2:end 1], :)./255, 'Label', [networks17 'Subcortex'])
pos = get(gcf, 'Position');
%set(gcf, 'Position', [0 0 pos(3:4)*2])
export_fig([outname '_neg'], '-png', '-nofontswap', '-a1');
set(gcf, 'color', 'w')
hgexport(f, [outname '_neg'])
close

end