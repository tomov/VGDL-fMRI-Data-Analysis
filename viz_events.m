% visualize events for a subject

% https://www.mathworks.com/help/database/ug/mongo.html#d117e86584
conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')

% https://www.mathworks.com/help/database/ug/mongo.find.html
runs = find(conn, 'runs', 'query', '{"subj_id": "26"}', 'limit', 10)


figure;

for r = 1:length(runs)
    run = runs(r)
    blocks = run.blocks;

    dur = round(run.end_time - run.scan_start_time + 5) * 1000;
    st = zeros(1,dur);
    en = zeros(1,dur);

    offs = 1;
    st(round(offs * 1000)) = 16;
    en(round((run.end_time - run.scan_start_time + offs) * 1000)) = -16;

    x = -offs*1000:(dur - 1000*offs - 1);

    for b = 1:length(blocks)
        block = blocks(b);
        instances = block.instances;
        
        st(round((block.start_time - run.scan_start_time + offs) * 1000)) = 8;
        en(round((block.end_time - run.scan_start_time + offs) * 1000)) = -8;

        for i = 1:length(instances)
            instance = instances(i);

            st(round((instance.start_time - run.scan_start_time + offs) * 1000)) = 4;
            en(round((instance.end_time - run.scan_start_time + offs) * 1000)) = -4;

            q = sprintf('{"subj_id": "26", "run_id": %d, "block_id": %d, "instance_id": %d}', r, b, i);
            plays = find(conn, 'plays', 'query', q, 'limit', 10);
        end
    end

    subplot(length(runs), 1, r);
    plot(x, st, 'color', 'green');
    hold on;
    plot(x, en, 'color', 'red');

end

