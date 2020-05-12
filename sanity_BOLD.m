glmodel = 3;
subj = 1;
mask = 'masks/sphere_glm3_theory_change_flag_48_10_32_r=4mm.nii';

EXPT = vgdl_expt();

% plot timecourses with vs. without whitening and filter
% from ccnl_get_activations and ccnl_get_residuals
%
[mask_format, mask, Vmask] = get_mask_format_helper(mask);
modeldir = fullfile(EXPT.modeldir,['model',num2str(rsa.glmodel)],['subj',num2str(subj)]);
load(fullfile(modeldir,'SPM.mat'));

Y = spm_data_read(SPM.xY.VY, find(mask));
KWY = spm_filter(SPM.xX.K,SPM.xX.W*Y);
KY = spm_filter(SPM.xX.K, Y);
WY = SPM.xX.W*Y;

Y = mean(Y,2);
KWY = mean(KWY,2);
KY = mean(KY,2);
WY = mean(WY,2);

figure;

subplot(4,1,1);
plot(Y);
legend({'raw BOLD'});
xlabel('TR');

subplot(4,1,2);
plot(KWY);
legend({'filtered & whitened BOLD'});
xlabel('TR');

subplot(4,1,3);
plot(KY);
legend({'filtered BOLD'});
xlabel('TR');

subplot(4,1,4);
plot(WY);
legend({'whitened BOLD'});
xlabel('TR');


% residuals

% from spm_spm.m
Xs        = spm_sp('Set',SPM.xX.X);    % X
Xs.X      = full(Xs.X);
assert(immse(SPM.xX.X, Xs.X) < 1e-10);

KWXs        = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X));    % KWX
KWXs.X      = full(KWXs.X);
assert(immse(SPM.xX.xKXs.X, KWXs.X) < 1e-10);

KXs        = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.X));    % KX
KXs.X      = full(KXs.X);

WXs        = spm_sp('Set',SPM.xX.W*SPM.xX.X);    % KWX
WXs.X      = full(WXs.X);

res = spm_sp('r',Xs,Y);
kwres = spm_sp('r',KWXs,KWY);
kres = spm_sp('r',KXs,KY);
wres = spm_sp('r',WXs,WY);

figure;

subplot(4,1,1);
plot(res);
legend({'raw residuals'});
xlabel('TR');

subplot(4,1,2);
plot(kwres);
legend({'filtered & whitened residuals'});
xlabel('TR');

subplot(4,1,3);
plot(kres);
legend({'filtered residuals'});
xlabel('TR');

subplot(4,1,4);
plot(wres);
legend({'whitened residuals'});
xlabel('TR');

% repeat but with official methods

figure;
i = 1;
for whiten = 0:1
    for filter = 0:1
        subplot(4,1,i);
        i = i+1;
        res = ccnl_get_residuals(EXPT, glmodel, mask, subj, whiten, filter);
        res = mean(res{1}, 2);
        plot(res);
    end
end

figure;
i = 1;
for whiten = 0:1
    for filter = 0:1
        subplot(4,1,i);
        i = i+1;
        A = ccnl_get_activations(EXPT, glmodel, mask, subj, whiten, filter);
        A = mean(A{1},2);
        plot(A);
    end
end
