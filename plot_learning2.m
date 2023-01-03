% using data from plot_learning.m (generated with separate runs)

close all
clear all

% determined number of steps

max_human_steps = 1703; % average # human steps

% load models (one subjects)
%load(fullfile(get_mat_dir(),'plot_learning_subj12.mat'), 'results');
%load(fullfile(get_mat_dir(),'plot_learning_subj12_ablations2.mat'), 'results');
load(fullfile(get_mat_dir(),'plot_learning_subj12_ablations2.mat'), 'results');
empa = process_cumwins(results(2).cumwins, max_human_steps);
dqn = process_cumwins(results(3).cumwins, max_human_steps);
empa_no_intrinsic = process_cumwins(results(4).cumwins, max_human_steps);
empa_no_IW = process_cumwins(results(5).cumwins, max_human_steps);
empa_eps_greedy = process_cumwins(results(6).cumwins, max_human_steps);


% load humans (all subjects)
load(fullfile(get_mat_dir(),'plot_learning_humans1-11_data.mat'), 'data');
data_1_11 = data;
load(fullfile(get_mat_dir(),'plot_learning_humans12-32_data.mat'), 'data');
data_12_32 = data;
data_1_11=data_1_11(:,1:max_human_steps);
data_12_32=data_12_32(:,1:max_human_steps);
humans = [data_1_11; data_12_32];
humans = humans(2:end,:); % we didn't count steps properly for subject 1 due to a bug



figure('pos', [370 477 500 600]) ;


subplot(2,1,1);

cmap = [0.4460 0.6740 0.1880;
        0 0.4470 0.7410;
        0.8500 0.3250 0.0980;
         0.7085    0.6669    0.8734;
        0.5809    0.1964    0.5266;
            0.9184    0.7308    0.1890;
        ];


hold on;
linewidth=1.5;

h0 = plot(humans', 'color', [0.9 0.9 0.9]);

m = mean(humans,1);
se = std(humans,1)/sqrt(size(humans,1));
h1 = plot(m, 'color', cmap(1,:),'linewidth',linewidth);
%hf = fill([1:length(m) flip(1:length(m))], [m+se flip(m-se)], cmap(1,:));
%set(hf,'facealpha',0.3,'edgecolor','none');

h2 = plot(empa','color',cmap(2,:),'linewidth',linewidth);
h3 = plot(dqn','color',cmap(3,:),'linewidth',linewidth);
h4 = plot(empa_no_intrinsic','color',cmap(4,:),'linewidth',linewidth);
h5 = plot(empa_no_IW','color',cmap(5,:),'linewidth',linewidth);
h6 = plot(empa_eps_greedy','color',cmap(6,:),'linewidth',linewidth);

hold off;

l=legend([h0(1);h1(1);h2(1);h3(1);h4(1);h5(1);h6(1)],...
    {'Human (individual)','Human (mean)','EMPA','DDQN',...
     'EMPA, no IR', 'EMPA, no IW', 'EMPA, \epsilon-greedy'});%,'Interpreter','latex');
%set(l,'position',[0.1333 0.7596 0.4150 0.1637]);
%set(l,'position',[0.1997    0.0476    0.4950    0.2132]);
set(l,'position',[0.1430 0.8051 0.3332 0.1873]);
xlabel('steps taken by agent');
ylabel('episodes won');
%title('Human and model learning');
xlim([0 max_human_steps]);


%title('Human and model learning');


%% statistics

t = 1:length(m);

% null distribution
human_slopes = [];
for s=1:size(humans,1)
    mdl=fitlm(t',humans(s,:)', 'intercept',false);
    slope=mdl.Coefficients.Estimate;
    human_slopes=[human_slopes,slope];
end
human_slopes

model_data = {empa, dqn, empa_no_intrinsic, empa_no_IW, empa_eps_greedy};
model_names = {'EMPA','DDQN','EMPA, no intrinsic rewards', 'EMPA, no iterative width', 'EMPA, eps-greedy'};
model_slopes = [];
for i=1:length(model_data)
    mdl=fitlm(t',model_data{i}', 'intercept',false);
    slope = mdl.Coefficients.Estimate
    %[p,h,stats] = signrank(human_slopes-empa_slope)
    %[h,p,ci,stats] = ttest(human_slopes-empa_slope)
    [h,p,ci,stats] = ttest2(human_slopes,slope);
    fprintf('Human vs. %s slopes: t(%d)=%.3f, p=%.4f, two-sample t-test\n',model_names{i},stats.df,stats.tstat,p);
    model_slopes(i)=slope;
end



%figure('pos', [570 477 400 200]) ;
subplot(4,1,3);

hold on;
histogram(human_slopes,'facecolor',cmap(1,:));

for i=1:length(model_data)
    plot([model_slopes(i) model_slopes(i)], [0 18],'color',cmap(i+1,:),'linewidth',linewidth);
end

%legend({'Human','EMPA','DQN'});
xlabel('slope');
ylabel('count');











