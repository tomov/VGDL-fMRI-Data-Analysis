% using data from plot_learning.m (generated with separate runs)

close all
clear all

% determined number of steps

max_human_steps = 1703; % average # human steps

% load models (one subjects)
load(fullfile(get_mat_dir(),'plot_learning_subj12.mat'), 'results');
empa = process_cumwins(results(2).cumwins, max_human_steps);
dqn = process_cumwins(results(3).cumwins, max_human_steps);


% load humans (all subjects)
load(fullfile(get_mat_dir(),'plot_learning_humans1-11_data.mat'), 'data');
data_1_11 = data;
load(fullfile(get_mat_dir(),'plot_learning_humans12-32_data.mat'), 'data');
data_12_32 = data;
data_1_11=data_1_11(:,1:max_human_steps);
data_12_32=data_12_32(:,1:max_human_steps);
humans = [data_1_11; data_12_32];
humans = humans(2:end,:); % we didn't count steps properly for subject 1 due to a bug

figure('pos', [370 477 400 400]) ;

cmap = [0.4460 0.6740 0.1880;
        0 0.4470 0.7410;
        0.8500 0.3250 0.0980];


hold on;

h0 = plot(humans', 'color', [0.7 0.7 0.7]);

m = mean(humans,1);
se = std(humans,1)/sqrt(size(humans,1));
h1 = plot(m, 'color', cmap(1,:),'linewidth',3);
%hf = fill([1:length(m) flip(1:length(m))], [m+se flip(m-se)], cmap(1,:));
%set(hf,'facealpha',0.3,'edgecolor','none');

h2 = plot(empa','color',cmap(2,:),'linewidth',3);
h3 = plot(dqn','color',cmap(3,:),'linewidth',3);

hold off;

l=legend([h0(1);h1(1);h2(1);h3(1)],{'Human (individual)','Human (mean)','EMPA','DDQN'});
set(l,'position',[0.1333 0.7596 0.4150 0.1637]);
xlabel('Steps taken by agent');
ylabel('Episodes won');
title('Human and model learning');
xlim([0 max_human_steps]);

