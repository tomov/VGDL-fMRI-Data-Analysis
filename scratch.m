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

subj_id = 181

query = sprintf('{"subj_id": "%d"}', subj_id)

subj = find(conn, 'subjects', 'query', query, 'limit', 10)

runs = find(conn, 'runs', 'query', '{"subj_id": "177"}', 'limit', 10)

plays = find(conn, 'plays', 'query', '{"subj_id": "177"}', 'limit', 10)

regressors = find(conn, 'regressors', 'query', '{"subj_id": "177"}', 'limit', 10)
