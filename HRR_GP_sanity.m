%load('mat/HRR_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=10_sigma_w=1.mat')
load('mat/HRR_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=10_sigma_w=1.000_norm=1.mat');

figure;

subplot(2,4,1);
imagesc(theory_kernel);
title('theory_kernel', 'interpreter', 'none');
colorbar;

subplot(2,4,2);
imagesc(sprite_kernel);
title('sprite_kernel', 'interpreter', 'none');
colorbar;

subplot(2,4,3);
imagesc(interaction_kernel);
title('interaction_kernel', 'interpreter', 'none');
colorbar;

subplot(2,4,4);
imagesc(termination_kernel);
title('termination_kernel', 'interpreter', 'none');
colorbar;


subplot(2,4,5);
imagesc(theory_Xx);
title('theory_Xx', 'interpreter', 'none');
colorbar;

subplot(2,4,6);
imagesc(sprite_Xx);
title('sprite_Xx', 'interpreter', 'none');
colorbar;

subplot(2,4,7);
imagesc(interaction_Xx);
title('interaction_Xx', 'interpreter', 'none');
colorbar;

subplot(2,4,8);
imagesc(termination_Xx);
title('termination_Xx', 'interpreter', 'none');
colorbar;

