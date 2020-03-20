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

subj_id = 1;
run_id = 5;
game_name = 'vgfmri3_chase'

limit = 100;

query = sprintf('{"subj_id": "%d"}', subj_id)
subj = find(conn, 'subjects', 'query', query, 'limit', limit)

query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id)
runs = find(conn, 'runs', 'query', query, 'limit', limit)

query = sprintf('{"subj_id": "%d", "run_id": %d, "game_name": "%s"}', subj_id, run_id, game_name)
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
