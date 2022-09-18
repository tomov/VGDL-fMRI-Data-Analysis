
% In response to neuron reviewer 1, comment about learning

function show_figure(figure_name)

switch figure_name

    case 'corr_theory_HRR_theory_update__subj_1'
        load(fullfile(get_mat_dir(false), 'corr_theory_HRR_theory_update.mat'));

        EXPT = vgdl_expt;
        glmodel = 102;
        subj_id = 1;
        [~, theory_update] = load_GLM_kernel(EXPT, glmodel, subj_id, {'theory_change_flag'}, false, false);
        load(sprintf('/n/holystore01/LABS/gershman_lab/Users/mtomov13/VGDL/mat/HRR_cannon_repro_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1_concat=0_novelty=1.mat', subj_id), 'theory_Xx', 'r_id');
        [coeff,score,latent,tsquared,explained,mu] = pca(theory_Xx);
        
        figure('position', [557 387 386 489]);
        
        subplot(3,1,1)
        hold on;
        plot(theory_Xx(283:2*283,1) * 10);
        plot(theory_update(283:2*283));
        legend({'theory (top PC)', 'theory update'});
        xlabel('time within run (TR)');
        ylabel('signal (a.u.)');
        title('Subject 1');

        subplot(3,1,2)
        bar(cumsum(mean(explaineds(1,1:30), 1)));
        xlabel('Top theory PCs');
        ylabel({'cumulative variance', 'explained (%)'});

        subplot(3,2,5)
        bar(rs(1,1:7));
        hold on;
        xlabel('Top theory PCs');
        ylabel('Pearson''s r');
        %title('                   corr(theory HRR PC, theory update)');

        subplot(3,2,6)
        bar(ps(1,1:7));
        hold on;
        xlabel('Top theory PCs');
        ylabel('p-value');
        %bar(ps(1,:));

        print('svg/neuron_revision/figure_neuron_R2_corr_theory_HRR_theory_update__subj_1.svg', '-dsvg');
        
    
    case 'corr_theory_HRR_theory_update__all_subj'
        load(fullfile(get_mat_dir(false), 'corr_theory_HRR_theory_update.mat'));
        
        zs = atanh(rs);
        zs = zs(:,1:30);
        [sem, me] = wse(zs);
        [h,p,ci,stats] = ttest(zs);
        p

        figure('position', [557 387 386 489]);

        subplot(3,1,1);
        bar(cumsum(mean(explaineds(:,1:30), 1)));
        hold on;
        xlabel('Top theory PCs');
        ylabel({'cumulative variance', 'explained (%)'});
        title('All subjects');

        subplot(3,1,2);
        hold on;
        bar(me);
        h = errorbar(me, sem, '.', 'MarkerSize', 1);
        h.CapSize = 1;
        %ax = gca;
        %for j = 1:size(zs, 2)
        %    if me(j) < 0
        %        y = me(j) - sem(j) - 0.01;
        %    else
        %        y = me(j) + sem(j) + 0.01;
        %    end
        %    h = text(j, y, sprintf('p = %.2f', p(j)), 'fontsize', 6, 'HorizontalAlignment', 'center');
        %    set(h, 'rotation', 60);
        %end
        xlabel('Top theory PCs');
        ylabel('z');
        ylim([-0.08 0.06]);

        subplot(3,1,3);
        bar(p);
        hold on;
        xlabel('Top theory PCs');
        ylabel('p-value');
        %title('                   corr(theory HRR PC, theory update)');

        print('svg/neuron_revision/figure_neuron_R2_corr_theory_HRR_theory_update__all_subj.svg', '-dsvg');


    otherwise
        assert(false, 'Invalid figure name');
end
