rsa_idx = 1;
subj = 1;
%mask = 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii';

load('mat/neurosynth_rsa_1_us=0_l=1_nperms=1000_nroi=351.mat');
mask = roi_masks{find(region == 307)};

EXPT = vgdl_expt_nosmooth();
rsa = vgdl_create_rsa(rsa_idx, subj);

figure;
plot(rsa.model(1).features);
title('features');
xlabel('time');

Behavioral = ccnl_behavioral_rdms(EXPT, rsa_idx, subj); % for plotting

figure;
imagesc(Behavioral(1).subj(subj).RDM);
colorbar;
title('model RDM');
xlabel('time');
ylabel('time');

    
rng('shuffle');
figure;
create_rsa = EXPT.create_rsa;
for i = 1:5
    seed = randi(1000000);
    EXPT.create_rsa = @(rsa_idx, subj_id) create_rsa(rsa_idx, subj_id, seed); % overwrite create_rsa with one that randomly shuffles 

    Behavioral_perm = ccnl_behavioral_rdms(EXPT, rsa_idx, subj);
    subplot(1,5,i);
    imagesc(Behavioral_perm(1).subj(subj).RDM);
end
colorbar;
title('permuted model RDM');


%EXPT = vgdl_expt_nosmooth();
Neural = ccnl_roi_rdms(EXPT, rsa_idx, mask, subj, true);

figure;
imagesc(Neural(1).subj(1).B);
xlabel('voxels');
ylabel('time');
colorbar;
title('neural activity');

figure;
imagesc(Neural(1).subj(1).RDM);
colorbar;
title('Neural RDM');
xlabel('time');
ylabel('time');



