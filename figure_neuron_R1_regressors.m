function show_figure(figure_name)

switch figure_name
    
    case 'theory_HRRs'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_Xx');

        figure('pos', [712 152 479 764]);
        imagesc(theory_Xx);
        title('EMPA theory HRRs: design matrix');
        xlabel('regressor');
        ylabel('TR');
        colorbar

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_HRRs.svg', '-dsvg');



    case 'theory_kernel'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_kernel');

        figure('pos', [712 152 764 764]);
        imagesc(theory_kernel);
        title('EMPA theory HRRs: GP kernel');
        xlabel('TR');
        ylabel('TR');
        colorbar

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_kernel.svg', '-dsvg');




    case 'theory_kernel_GLM1_projected'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_kernel');

        EXPT = vgdl_expt;
        glmodel = 1;
        subj = 1;
        %[mask_format, mask, Vmask] = get_mask_format_helper('masks/mask.nii');
        %[Y, K, W, R, SPM_run_id, R_] = load_BOLD(EXPT, glmodel, subj, mask, Vmask);

        % TODO dedupe with load_BOLD.m
        % get K, W, R
        modeldir = fullfile(EXPT.modeldir,['model',num2str(glmodel)],['subj',num2str(subj)]);
        load(fullfile(modeldir,'SPM.mat'));
        X = SPM.xX.X; % original design matrix
        K = SPM.xX.K; % high-pass filter
        W = SPM.xX.W; % whitening matrix
        KWX = SPM.xX.xKXs.X; % high-pass filtered & whitened design matrix

        R = spm_sp('r',SPM.xX.xKXs); % residual forming matrix R = I - X * pinv(X)

        % convert filter K to matrix form
        % see spm_filter.m
        for s = 1:length(K)
            I(K(s).row,K(s).row) = eye(length(K(s).row));
            X0X0(K(s).row, K(s).row) = K(s).X0*K(s).X0';
        end
        K = I - X0X0; % high-pass filter matrix

        assert(immse(K*W*X, KWX) < 1e-15);
        assert(immse(R, eye(size(X,1)) - SPM.xX.xKXs.u*SPM.xX.xKXs.u') < 1e-15);


        ker = R*K*W*theory_kernel*W'*K'*R';

        figure('pos', [712 152 764 764]);
        imagesc(ker);
        title('EMPA theory HRRs: GP kernel, controlling for game identity');
        xlabel('TR');
        ylabel('TR');

        set(gca,'ColorScale','log');
        h = colorbar;
        set(h, 'XTickLabel', arrayfun(@(x) sprintf('%.4f', x), exp(h.Ticks), 'UniformOutput', false));

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_kernel_GLM1_projected.svg', '-dsvg');



    otherwise
        assert(false, 'Invalid figure name');
end
