clear

%% extract filter matrix from SPM

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
R1 = spm_sp('r',SPM.xX.xKXs); % residual forming matrix R = I - X X^+

u = SPM.xX.xKXs.u;
R2 = eye(size(X,1)) - u*u';
assert(immse(R1, R2) < 1e-20);

%





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
