function ABCD_LRR_collect_opt_param(model_dir, bhvr_ls)

% ABCD_LRR_collect_opt_param(model_dir, bhvr_ls)
%
% 

ls_dir = fullfile(getenv('HOME'), 'storage', 'MyProject', 'fairAI', 'ABCD_race', 'scripts', 'lists');
if(~exist('bhvr_ls', 'var') || isempty(bhvr_ls))
    bhvr_ls = fullfile(ls_dir, 'behavior_list.txt');
end

[bhvr_nm, nbhvr] = CBIG_text2cell(bhvr_ls);
Nfolds = 120;

opt_lambda = zeros(nbhvr, Nfolds);
for b = 1:nbhvr
    for f = 1:Nfolds
        opt = load(fullfile(model_dir, bhvr_nm{b}, 'params', ['fold_' num2str(f)], ...
            ['selected_parameters_' bhvr_nm{b} '.mat']));
        opt_lambda(b,f) = opt.curr_lambda;
    end
end

mode_lambda = mode(opt_lambda, 2);
med_lambda = median(opt_lambda, 2);

save(fullfile(model_dir, 'opt_lambda.mat'), 'opt_lambda', 'mode_lambda', 'med_lambda')

    
end