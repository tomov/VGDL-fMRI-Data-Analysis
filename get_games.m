function [game_names, onsets, durs] = get_games(subj_id, run, conn)

    % helper function to get game/instance boxcar regressors in vgdl_create_multi
    % copied & modified from GLM 1
    % note run is a struct
    %

    game_names = {};
    onsets = {};
    durs = {};

    blocks = run.blocks;
    for b = 1:length(blocks)
        block = blocks(b);
        instances = block.instances;

        game_name = block.game.name;
        game_names = [game_names; {game_name}];
        
        ons = [];
        dur = [];
        for i = 1:length(instances)
            instance = instances(i);

            st = instance.start_time - run.scan_start_ts;
            en = instance.end_time - run.scan_start_ts; % includes "YOU WON / LOST / etc" screen

            ons = [ons, st];
            dur = [dur, en - st];
        end

        onsets = [onsets; {ons}];
        durs = [durs; {dur}];
    end
