function [Pset, B] = ABCD_multilevel_block_perm(site, family, grp, Nperm)

    % Pset = ABCD_multilevel_block_perm(site, family, rel_grp, Nperm)
    %
    % 1. Generate multi-level block definitions B for palm_tree package.
    % 2. Call palm_tree to create graph representation of exchangeable blocks.
    % 3. Call palm_permtree to generate permuted indices.
    %
    % First level: whole dataset (-/+: whether sites can be permuted, always negative)
    % Second level: sites (-/+: whether types of families can be permuted, always negative)
    % Third level: type of families, i.e. families with the same structures, 
    %               e.g. same number of non-twin siblings, same number of twins with 2 people, ...
    %              (-/+: whether families of the same type can be permuted, always positive)
    % Fourth level: families (-/+: whether groups can be permuted)
    % Fifth level: twins / within-family groups 
    %              (-/+: whether individuals in the same group can be permuted, always positive)
    % Sixth level: individuals (always positive)
    
    addpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', 'non_default_packages', 'palm', 'palm-alpha109'))
    
    %% generate block definations
    B = -ones(length(site), 1);      % first level: unexchangeable across sites
    all_idx = [1:length(site)]';
    
    uniq_sites = unique(site);
    count_site = 1;
    for s = 1:length(uniq_sites)
        idx_currsite = all_idx(strcmp(site, uniq_sites(s)));
        B(idx_currsite, 2) = -count_site;   % each site has a block number at the second level
    
        curr_fam = family(idx_currsite);
        uniq_fam = unique(curr_fam);
    
        N_perfam = zeros(length(uniq_fam), 1);
        for f = 1:length(uniq_fam)
            N_perfam(f) = length(find(curr_fam == uniq_fam(f)));
        end
    
        % column1: scalar, number of non-twin siblings in this family; 
        % column2: scalar, number of twin groups with 2 subjects
        % column3: scalar, number of twin groups with 3 subjects; ...
        fam_struct = zeros(length(uniq_fam), max(N_perfam));
        has_singleton = 0;
        for f = 1:length(uniq_fam)
            idx_currfam = idx_currsite(curr_fam == uniq_fam(f));
            curr_grp = grp(idx_currfam);
            uniq_grp = unique(curr_grp);
            N_pergrp = zeros(length(uniq_grp), 1);
            for g = 1:length(uniq_grp)
                N_pergrp(g) = length(find(curr_grp == uniq_grp(g)));
            end
    
            for n = 1:max(N_perfam)
                fam_struct(f, n) = length(find(N_pergrp == n));
            end
            
            if(~any(fam_struct(f, 2:end)) && fam_struct(f,1) == 1)    % singleton
                B(idx_currfam, 3:5) = 1;        % singletons are defined as block 1 in the third to the fifth levels
                has_singleton = 1;
            end
        end
    
        [uniq_famstr, ia, ic] = unique(fam_struct, 'rows');
        if(has_singleton)
            count_str = 2;   % all non-singleton families are numbered starting from 2
        else
            count_str = 1;
        end
        for r = 1:size(uniq_famstr, 1)
            if(~any(uniq_famstr(r, 2:end)) && uniq_famstr(r, 1) == 1)
                continue   % skip all singletons because they are assigned in the for-loop above
            end
    
            if(length(find(uniq_famstr(r,:))) == 1)
                % if this type of families only have 1 type of childrens 
                % (e.g. only non-twins, or only 2-people twins, or only triplets, ...), 
                % groups can be permuted
                signed = 1;
            else
                % if there are multiple types of groups in this type of families, 
                % these groups cannot be permuted at 4th level (family level), 
                % subjects can only be permuted within twins (5th level)
                signed = -1;
            end
    
    
            fam_idx = find(ic==r);   % row indices of the families with the same structure as the current loop's structure
            count_fam = 1;
            for i = 1:length(fam_idx)
                idx_hold = idx_currsite(curr_fam == uniq_fam(fam_idx(i)));
                B(idx_hold, 3) = count_str;
    
                B(idx_hold, 4) = signed * count_fam;
    
                grp_hold = grp(idx_hold);
                uq_grp_hold = unique(grp_hold);
                count_grp = 1;
                for j = 1:length(uq_grp_hold)
                    B(idx_hold(grp_hold == uq_grp_hold(j)), 5) = count_grp;
                    count_grp = count_grp + 1;
                end
    
                count_fam = count_fam + 1;
            end
            count_str = count_str + 1;
        end
        count_site = count_site + 1;
    end
    B(:, 6) = 1:length(site);
    
    %% call PALM package
    Ptree = palm_tree(B);
    Pset = palm_permtree(Ptree, Nperm, false, true);
    
    rmpath(fullfile(getenv('CBIG_CODE_DIR'), 'external_packages', 'matlab', 'non_default_packages', 'palm', 'palm-alpha109'))
        
    end