clear

% https://www.mathworks.com/help/database/ug/mongo.html#d117e86584
conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')

% https://www.mathworks.com/help/database/ug/mongo.find.html
%{
%docs = find(conn, 'states')
docs = find(conn, 'states', 'query', '{"exp_id": "ry44FBXnr", "game_name": "vgfmri2_helper", "game_round": "1"}', 'limit', 10)

%val = jsondecode(docs(10).game_states)
val = jsondecode(docs(10).game_real_states);

val
%}

%docs = find(conn, 'subjects', 'query', '{"subj_id": "21"}', 'limit', 10)

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
