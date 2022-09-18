clear all;

EXPT = vgdl_expt;
glmodel = 102;

[allSubjects, subj_dirs, goodRuns, goodSubjs] = vgdl_getSubjectsDirsAndRuns();

explaineds = [];
rs = [];
ps = [];

for subj_id = 1:32
    subj_id

    run_ids = find(goodRuns{find(subj_id == allSubjects)});

    % Theory update
    [~, theory_update] = load_GLM_kernel(EXPT, glmodel, subj_id, {'theory_change_flag'}, false, false);

    % Theory representation
    load(sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', subj_id), 'theory_Xx', 'r_id');
    which_TRs = ismember(r_id, run_ids);
    theory_Xx = theory_Xx(which_TRs, :);


    [coeff,score,latent,tsquared,explained,mu] = pca(theory_Xx);
    explaineds = [explaineds; explained'];

    [r, p] = corr(score(:,1:100), theory_update);
    rs = [rs; r'];
    ps = [ps; p'];
end

save(fullfile(get_mat_dir(false), 'corr_theory_HRR_theory_update.mat'));
