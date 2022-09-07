function run = get_run(subj_id, run_id)

    conn = mongo('holy7c22108.rc.fas.harvard.edu', 27017, 'heroku_7lzprs54', 'UserName', 'reader', 'Password', 'parolatamadafaqa');

    query = sprintf('{"subj_id": "%d", "run_id": %d}', subj_id, run_id); % in python we index runs from 0 (but not subjects) 

    run = find(conn, 'runs', 'query', query);
    assert(length(run) == 1);
