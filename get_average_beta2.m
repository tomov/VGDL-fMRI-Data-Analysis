% get average beta from GLM across a bunch of subjects
% separate runs

close all;
clear all;

EXPT = vgdl_expt();


train_subjects = 1:2:length(EXPT.subject);

[whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');

load(fullfile(get_mat_dir(false), 'glm_bic_bms_atlas=AAL2_GLM_102_multiplex_with_controls.mat'), 'glms');
glmodels = glms;

for glmodel = glmodels

    fprintf('  glm %d\n', glmodel);
    
    tic
    disp('           loading initial betas...');
    [B, names] = ccnl_get_beta_series(EXPT, glmodel, 1, ' ', whole_brain_mask); % all
    toc

    cum_B = [];

    for s = 1:length(train_subjects)
        subj_id = train_subjects(s);

        tic
        fprintf('           loading betas for %d...\n', subj_id);
        [B, names] = ccnl_get_beta_series(EXPT, glmodel, subj_id, ' ', whole_brain_mask); % all
        toc

        if isempty(cum_B)
            cum_B = B;
            cnt = ones(size(B,1),1);
            avg_names = names;
        else
            % see if there's any new regressors
            ix = find(~ismember(names, avg_names));
            if length(ix) > 0
                cum_B = [cum_B; zeros(length(ix), size(B,2))];
                cnt = [cnt; zeros(length(ix), 1)];
                avg_names = [avg_names; names(ix)];
            end

            for i = 1:size(B,1)
                ix = find(strcmp(names{i}, avg_names));
                cum_B(ix,:) = cum_B(ix,:) + B(i,:);
                cnt(ix,:) = cnt(ix,:) + 1;
            end
        end

    end

    avg_B = cum_B ./ cnt;

    filename = sprintf('get_average_beta2_glm=%d.mat', glmodel);
    fprintf('Saving to %s\n', filename);
    save(fullfile(get_mat_dir(2), filename));
end
