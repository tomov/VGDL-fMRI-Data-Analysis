% visualize events for a subject

% https://www.mathworks.com/help/database/ug/mongo.html#d117e86584
conn = mongo('127.0.0.1', 27017, 'heroku_7lzprs54')

subj_id = 181;

% https://www.mathworks.com/help/database/ug/mongo.find.html
runs = find(conn, 'runs', 'query', sprintf('{"subj_id": "%d"}', subj_id), 'limit', 10)


figure;

for r = 1:length(runs)
    run = runs(r);
    blocks = run.blocks;

    subplot(length(runs), 1, r);
    hold on;

    dur = round(run.end_time - run.scan_start_ts + 5) * 1000;
    rr = zeros(1,dur); % runs
    bb = zeros(1,dur); % blocks
    ii = zeros(1,dur); % instances
    pp = zeros(1,dur); % plays
    ee = zeros(1,dur); % events
    aa = zeros(1,dur); % actions

    offs = 1;
    st = round(offs * 1000) ;
    en = round((run.end_time - run.scan_start_ts + offs) * 1000);
    rr(st:en) = 8;
    text((st*0.7+en*0.3), 7, sprintf('run %d', run.run_id), 'interpreter', 'none');

    x = -offs*1000:(dur - 1000*offs - 1);

    for b = 1:length(blocks)
        block = blocks(b);
        instances = block.instances;
       
        st = round((block.start_time - run.scan_start_ts + offs) * 1000) + 3;
        en = round((block.end_time - run.scan_start_ts + offs) * 1000) - 3;
        bb(st:en) = 6;
        text((st*0.7+en*0.3), 5, block.game.name, 'interpreter', 'none');

        for i = 1:length(instances)
            instance = instances(i);

            st = round((instance.start_time - run.scan_start_ts + offs) * 1000) + 5;
            en = round((instance.end_time - run.scan_start_ts + offs) * 1000) - 5;
            ii(st:en) = 4;
            text((st*0.7+en*0.3), 3, sprintf('level %d', instance.level_id), 'interpreter', 'none');

            q = sprintf('{"subj_id": "%d", "run_id": %d, "block_id": %d, "instance_id": %d}', subj_id, run.run_id, block.block_id, instance.instance_id);
            plays = find(conn, 'plays', 'query', q);

            for p = 1:length(plays)
                play = plays(p);

                st = round((play.start_time - run.scan_start_ts + offs) * 1000) + 7;
                en = round((play.end_time - run.scan_start_ts + offs) * 1000) - 7;
                pp(st:en) = 2;
                if isempty(play.win)
                    win = 'T';
                elseif play.win
                    win = 'W';
                else
                    win = 'L';
                end
                text((st*0.9+en*0.1), 1, sprintf('%s', win), 'interpreter', 'none');

                for a = 1:length(play.actions)
                    aa(round((play.actions{a}{2} - run.scan_start_ts + offs) * 1000)) = -0.2;
                end

                for e = 1:length(play.events)
                    ee(round((play.events(e).ts - run.scan_start_ts + offs) * 1000)) = 0.2;
                end
            end
        end
    end

    title(['Subject ', subj_id]);

    plot(x, rr);
    plot(x, bb);
    plot(x, ii);
    plot(x, pp);
    plot(x, aa);
    plot(x, ee);
    legend({'run', 'blocks', 'instances', 'plays', 'actions', 'events'});

end

