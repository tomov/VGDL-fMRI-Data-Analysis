function show_figure(figure_name)

figure_scale = 0.7;

switch figure_name
    
    case 'theory_HRRs'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_Xx');

        figure('pos', [712 152 479*figure_scale 764*figure_scale]);
        imagesc(theory_Xx);
        title('EMPA theory HRRs: design matrix');
        xlabel('regressor');
        ylabel('TR');
        colorbar

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_HRRs.svg', '-dsvg');



    case 'theory_kernel'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_kernel');

        figure('pos', [712 152 764*figure_scale 764*figure_scale]);
        imagesc(theory_kernel);
        title('EMPA theory HRRs: GP kernel');
        xlabel('TR');
        ylabel('TR');
        colorbar

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_kernel.svg', '-dsvg');

    case 'tsne'


    %
    % --------- same but GLM1 projected -------------
    %
    
    case 'theory_HRRs_GLM1_projected'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_Xx');

        EXPT = vgdl_expt;
        glmodel = 1;
        subj = 1;
        [R, K, W] = get_R_K_W(EXPT, glmodel, subj);


        figure('pos', [712 152 479*figure_scale 764*figure_scale]);
        imagesc(R*K*W*theory_Xx);
        title('EMPA theory HRRs: design matrix');
        xlabel('regressor');
        ylabel('TR');

        set(gca,'ColorScale','log');
        h = colorbar;
        %set(h, 'XTickLabel', arrayfun(@(x) sprintf('%.4f', x), exp(h.Ticks), 'UniformOutput', false));

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_HRRs.svg', '-dsvg');




    case 'theory_kernel_GLM1_projected'

        load('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=1_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', 'theory_kernel');

        EXPT = vgdl_expt;
        glmodel = 1;
        subj = 1;
        [R, K, W] = get_R_K_W(EXPT, glmodel, subj);
        ker = R*K*W*theory_kernel*W'*K'*R';

        figure('pos', [712 152 764*figure_scale 764*figure_scale]);
        imagesc(ker);
        title('EMPA theory HRRs: GP kernel, controlling for game identity');
        xlabel('TR');
        ylabel('TR');

        set(gca,'ColorScale','log');
        h = colorbar;
        %set(h, 'XTickLabel', arrayfun(@(x) sprintf('%.4f', x), exp(h.Ticks), 'UniformOutput', false));

        print('svg/neuron_revision/figure_neuron_R1_regressors_theory_kernel_GLM1_projected.svg', '-dsvg');




    otherwise
        assert(false, 'Invalid figure name');
end
