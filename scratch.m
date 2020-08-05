

%
%% sanity check to run after decode_gp_CV
%

load decode_gp_CV.mat

% load BOLD
use_smooth = true;
glmodel = 9;
maskfile = 'masks/ROI_x=42_y=28_z=26_1voxels_Sphere1.nii';
[mask_format, mask, Vmask] = get_mask_format_helper(maskfile);
[Y, K, W, R, run_id_TRs] = load_BOLD(EXPT, glmodel, subj_id, mask, Vmask);
Y = R*K*W*Y;
assert(size(Y,2) == 1); % single voxel
y = Y;

% fit GP w/ CV to get r
[theory_kernel, ~, HRRs, Xx] = gen_kernel_from_theory_id_seq(unique_theory_HRRs, theory_id_seq_orig, ts, run_id_frames);
ker = R*K*W*theory_kernel*W'*K'*R';
[r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, Y, ker, run_id_TRs);

% fit manually too to get LME
hyp = hyp_from_ker(nearestSPD(ker));
hyp.lik = log(1); % TODO sigma_n = 1 const
nlz = gp(hyp, @infGaussLik, meanfun, covfun, likfun, x, y);
logmarglik = -nlz;


% compare with pre-computed 
%
r_CV_1 = r_CV;
R2_CV_1 = R2_CV;
logmarglik_1 = logmarglik;
load('mat/fit_gp_CV_HRR_subj=1_us=1_glm=9_mask=mask_theory.mat', 'r_CV', 'R2_CV', 'mask');
whole_brain_mask = mask;
[mask] = ccnl_load_mask(maskfile);
which = mask(whole_brain_mask);
r_CV = r_CV(:,which);
R2_CV = R2_CV(:,which);

r_CV_1
r_CV

logmarglik_1
logmarglik


%{



load(sprintf('mat/HRR_subject_kernel_subj=%d_K=10_N=10_E=0.050_nsamples=100_sigma_w=1.000_norm=1.mat', subj_id), 'theory_kernel');
[r_CV, R2_CV, MSE_CV, SMSE_CV] = fit_gp_CV_simple(subj_id, use_smooth, glmodel, maskfile, theory_kernel);

r_CV_2 = r_CV;
R2_CV_2 = R2_CV;


for subj_id = 1:8
    %r_CV_1 = r_CV;
    %R2_CV_1 = R2_CV;
    load(sprintf('mat/fit_gp_CV_HRR_subj=%d_us=1_glm=9_mask=mask_theory.mat', subj_id), 'r_CV', 'R2_CV', 'mask');
    whole_brain_mask = mask;
    [mask] = ccnl_load_mask(maskfile);
    which = mask(whole_brain_mask);
    r_CV = r_CV(:,which);
    R2_CV = R2_CV(:,which);

    rs(subj_id,:) = mean(r_CV,1);
end

% Fisher z transform
zs = atanh(rs);

% t-test Pearson corr across subjects
[h,p,ci,stats] = ttest(zs);
ts = stats.tstat;

%}


%{
%% get voxel coordinates for Daphne

[mask, Vmask] = ccnl_load_mask('masks/mask_nosmooth.nii');

[x y z] = ind2sub(size(mask), find(mask));

cor = [x y z];
mni = cor2mni(cor, Vmask.mat);

save('mat/coords_nosmooth.mat');
%}


%{
clear

%% extract K, W, and R matrices from SPM

% svd sanity
%
X = rand(30,10);
[u,s,v] = svd(X,0);
assert(immse(u*s*v', X) < 1e-20);


load('mat/SPM.mat');

X = SPM.xX.X;
K = SPM.xX.K;
W = SPM.xX.W;
KWX = SPM.xX.xKXs.X;

%}

% spm_filter.m
%
%{
KX1 = spm_filter(K, X);

clear I;
clear X0X0;
for s = 1:length(K)
    I(K(s).row,K(s).row) = eye(length(K(s).row));
    X0X0(K(s).row, K(s).row) = K(s).X0*K(s).X0';
end

K = I - X0X0;
KX2 = K * SPM.xX.X;
assert(immse(KX1, KX2) < 1e-20);
%}

% residual forming matrix
%
%{
R1 = spm_sp('r',SPM.xX.xKXs); % residual forming matrix R = I - X X^+

u = SPM.xX.xKXs.u;
R2 = eye(size(X,1)) - u*u';
assert(immse(R1, R2) < 1e-20);

%}





%% generate searchlight ROIs with different params
%{

r = [4 6 10] / 1.5;

for use_smooth = 0:1
    for i = 1:length(r)
        radius = r(i);

        if use_smooth
            EXPT = vgdl_expt();
            [whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask.nii');
        else
            EXPT = vgdl_expt_nosmooth();
            [whole_brain_mask, Vwhole_brain_mask] = ccnl_load_mask('masks/mask_nosmooth.nii');
        end

        ROI = get_searchlight_rois(whole_brain_mask, Vwhole_brain_mask, radius);

        save(sprintf('mat/get_searchlight_rois_us=%d_r=%.4f.mat', use_smooth, radius), '-v7.3');
    end
end
%}

% https://www.mathworks.com/help/database/ug/mongo.html#d117e86584
%conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')

% https://www.mathworks.com/help/database/ug/mongo.find.html
%{
%docs = find(conn, 'states')
docs = find(conn, 'states', 'query', '{"exp_id": "ry44FBXnr", "game_name": "vgfmri2_helper", "game_round": "1"}', 'limit', 10)

%val = jsondecode(docs(10).game_states)
val = jsondecode(docs(10).game_real_states);

val
%}

%docs = find(conn, 'subjects', 'query', '{"subj_id": "21"}', 'limit', 10)

%{
subj_id = 2;
run_id = 1;
%game_name = 'vgfmri3_chase'

limit = 100;

query = sprintf('{"subj_id": "%d"}', subj_id)
subj = find(conn, 'subjects', 'query', query, 'limit', limit)

query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id)
runs = find(conn, 'runs', 'query', query, 'limit', limit)

query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id)
plays = find(conn, 'plays', 'query', query, 'limit', limit)

regressors = find(conn, 'regressors', 'query', query, 'limit', limit)
%}

%for r = 1:length(subj.runs)
%    run = subj.runs(r);
%    fprintf('run %d\n', r);
%    for b = 1:length(run.blocks)
%        block = run.blocks(b);
%        fprintf('   %s\n', block.game.name);
%    end
%end
%



%% count plays

%tot = 0;
%for s = 1:1
%    query = sprintf('{"subj_id": "%d"}', s)
%    subj = find(conn, 'subjects', 'query', query);;;;
%    for r = 1:length(subj.runs)
%        run = subj.runs(r);
%        for b = 1:length(run.blocks)
%            block = run.blocks(b);
%            for i = 1:length(block.instances)
%                instance = block.instances(i);
%                query = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d}', s, run.run_id, block.block_id, instance.instance_id);
%                cnt = count(conn, 'plays', 'query', query);
%                fprintf('%s => %d\n', query, cnt);
%                tot = tot + cnt;
%            end
%        end
%    end
%
%    tot
%end
%
%tot
