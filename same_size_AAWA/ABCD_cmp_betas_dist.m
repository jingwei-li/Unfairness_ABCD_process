function ABCD_cmp_betas_dist(fig_dir, orig_dir, resample_dir, bhvr_ls, colloq_ls, cfds_ls, Nfolds)

% ABCD_cmp_betas_dist(fig_dir, orig_dir, resample_dir, bhvr_ls, colloq_ls, cfds_ls, Nfolds)
%
% 

ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end
if(~exist('colloq_ls', 'var') || isempty(colloq_ls))
    colloq_ls = fullfile(ls_dir, 'colloquial_list.txt');
end
if(~exist('cfds_ls', 'var') || isempty(cfds_ls))
    cfds_ls = fullfile(ls_dir, 'confounds_list.txt');
end
if(~exist('Nfolds', 'var') || isempty(Nfolds))
    Nfolds = 120;
end

[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
colloq_nm = CBIG_text2cell(colloq_ls);
[cfds_nm, ncfds] = CBIG_text2cell(cfds_ls);
idx = strcmp(cfds_nm, 'peduc_avg');
if(any(idx))
    cfds_nm{idx} = 'Parental education';
end
N_bstrp = 100;

for b = 1:nbhvr
    orig_betas = zeros(Nfolds, ncfds+1);
    resampled_betas = zeros(N_bstrp, Nfolds, ncfds+1);
    for f = 1:Nfolds
        orig = load(fullfile(orig_dir, bhvr_nm{b}, 'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']));
        orig_betas(f, :) = orig.beta;

        for i = 1:N_bstrp
            resampled = load(fullfile(resample_dir, ['sample' num2str(i)], bhvr_nm{b}, ...
                'y', ['fold_' num2str(f)], ['y_regress_' bhvr_nm{b} '.mat']));
            resampled_betas(i, f, :) = resampled.beta;
        end
    end

    outname = fullfile(fig_dir, bhvr_nm{b});
    vis_two_beta_dist(orig_betas, squeeze(mean(resampled_betas, 1)), outname, ...
        [{'Constant 1'} cfds_nm], sprintf('Beta (%s)', colloq_nm{b}));
end
    
end




function vis_two_beta_dist(data_orig, data_resample, outname, cfds_nm, y_label)

% data_orig: Nfolds x #confounds
% data_resample: Nfolds x #confounds
    
color_orig = [80 80 80]./255;
color_resample = [130 130 130]./255;
legends = {'Original', 'Resampled'};

f = figure;
set(gcf, 'position', [0 0 900 350])
for c = 1:length(cfds_nm)
    subplot(1, length(cfds_nm), c);
    h1 = distributionPlot(data_orig(:,c), 'widthDiv', [2 1], 'histOri', 'left', 'color', color_orig, 'showMM', 6);
    subplot(1, length(cfds_nm), c);
    h2 = distributionPlot(gca, data_resample(:,c), 'widthDiv', [2 2], 'histOri', 'right', 'color', color_resample, 'showMM', 6);
    set(h1{2}(:), 'color', [211 94 96]./255)  % change color for 25%, 50%, 75% lines
    set(h2{2}(:), 'color', [211 94 96]./255)

    set(gca, 'tickdir', 'out', 'box', 'off')
    set(gca, 'xticklabel', cfds_nm{c}, 'fontsize', 12, 'linewidth', 2);
    rotateXLabels( gca(), 45 );
end

l = legend([h1{1}, h2{1}], legends);
set(l, 'fontsize', 11, 'linewidth', 2, 'location', 'best', 'box', 'off')

export_fig(outname, '-png', '-nofontswap', '-a1');
set(gcf, 'color', 'w');
hgexport(f, outname)
close

end