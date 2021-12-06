close all;
clear all;

% generated with semopy_bic.py in py_vgdl
filenames = {
    %fullfile(get_mat_dir(false), 'semopy_GLM_109_lmes.mat'),
    %fullfile(get_mat_dir(false), 'semopy_GLM_120_lmes.mat'),
    '/n/home_fasse/mtomov13/py_vgdl/mat/semopy_GLM_109_lmes.mat',
    '/n/home_fasse/mtomov13/py_vgdl/mat/semopy_GLM_120_lmes.mat',
    '/n/home_fasse/mtomov13/py_vgdl/mat/semopy_GLM_121_lmes.mat',
    '/n/home_fasse/mtomov13/py_vgdl/mat/semopy_GLM_122_lmes.mat',
}

event_names = {'theory update', 'theory update + 2 s', '+ 4 s', '+ 6 s'}

SEM_names = {'top down', 'bottom up'};

figure;

for i = 1:length(filenames)
    filenames{i}
    event_names{i}

    load(filenames{i}, 'lmes', 'bics');

    [alpha, exp_r, xp, pxp, bor] = bms(lmes);
    pxps(i,:) = pxp;
    bors(i,:) = bor;

    [sem, me] = wse(bics);

    subplot(1,length(filenames),i);
    hold on;
    bar(me - me(1));
    errorbar(me - me(1), sem, 'o', 'MarkerSize', 1);

    xlabel('SEM');
    ylabel('\Delta BIC');
    %xticklabels(glm_names(glm_ix));
    xticklabels(SEM_names);
    xtickangle(30);
    xticks(1:length(SEM_names));
    title(event_names{i}, 'interpreter', 'none');
end


table(event_names', pxps, bors)



figure;

ms = [];
sems = [];

for i = 1:length(filenames)
    load(filenames{i}, 'lmes', 'bics');
    dbic = bics(:,2) - bics(:,1);
    %[sems(i), ms(i)] = wse(dbic');
    ms(i) = mean(dbic);
    sems(i) = std(dbic) / sqrt(length(dbic));
end

hold on;
bar(ms);
errorbar(ms, sems, 'o', 'MarkerSize', 1);

xlabel('theory update offset');
ylabel('BIC_{bottom up} - BIC_{top down}');
%xticklabels(glm_names(glm_ix));
xticklabels({'+ 0 s', '+ 2 s', '+ 4 s', '+ 6 s'});
%xtickangle(30);
xticks(1:4);
title('SEM model comparison');
