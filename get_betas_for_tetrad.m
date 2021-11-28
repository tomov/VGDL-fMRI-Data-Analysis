% Extract Beta series for functional connectivity analysis (using TETRAD / IMaGES)

close all;
clear all;

EXPT = vgdl_expt();

atlas = 'AAL2_TETRAD_GLM_109';

fasse_ncf = false;
filename = fullfile(get_mat_dir(fasse_ncf), sprintf('get_betas_for_tetrad_atlas=%s.mat', atlas));
filename

%% get masks
%
[mask_filenames, roi_names, roi_masks] = get_anatomical_masks(atlas);
nROIs = length(mask_filenames)

subjects = 1:length(EXPT.subject);
nsubjects = length(subjects)

for s = 1:nsubjects
    s
    subj_id = subjects(s);
        
    for m = 1:nROIs
        B_tcf{s}(:, m) = mean(ccnl_get_beta_series(EXPT, 109, subj_id, 'theory_change_flag', mask_filenames{m}), 2);
    end
    
    % convert to tables
    B_tcf{s} = array2table(B_tcf{s}, 'VariableNames', roi_names);
end

filename
save(filename, '-v7.3');


for s = 1:nsubjects
    subj_id = subjects(s);
    table_filename = fullfile(get_mat_dir(fasse_ncf), sprintf('get_betas_for_tetrad_theory_change_flag_SPMsubj=%d.txt', subj_id));
    table_filename
    writetable(B_tcf{s}, table_filename, 'Delimiter', '\t');
end
