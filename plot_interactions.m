clear all;
close all;

%[all_game_names, sprite_types_ordered] = get_game_names_ordered(12);
%subj_ids = 12:32;
[all_game_names, sprite_types_ordered] = get_game_names_ordered(11);
subj_ids = 1:11;

count_walls = true;

run_ids = 1:6;
all_levels = 1:9;


agents(1).name = 'Human';
agents(1).tag = '';
agents(2).name = 'Random';
agents(2).tag = 'attempt_1_states';
%agents(3).name = 'EMPA';
%agents(3).tag = 'attempt_1_states';

%% get data


for g = 1:length(all_game_names)
    game_name = all_game_names{g};
    game_name
    sprite_types = sprite_types_ordered{g};
        
    for a = 1:length(agents)
        agent_name = agents(a).name;
        agent_tag = agents(a).tag;
        agent_name

        for s = 1:length(subj_ids)
            subj_id = subj_ids(s);

            if strcmp(agent_name, 'Human')
                %filename = sprintf('mat/fmri_avatar_interactions_subj=%d.mat', subj_id);  % <-- interactions according to the rules
                filename = sprintf('mat/fmri_avatar_collisions_subj=%d.mat', subj_id); % <-- all collisions
            else
                %filename = sprintf('mat/fmri_agent_avatar_interactions_subj=%d_agent=%s_tag=%s.mat', subj_id, agent_name, agent_tag);
                filename = sprintf('mat/fmri_agent_avatar_collisions_subj=%d_agent=%s_tag=%s.mat', subj_id, agent_name, agent_tag);
            end
            load(filename, 'levels', 'game_names', 'object_types');
            object_types = cellstr(object_types);
            game_names = cellstr(game_names);
            levels = levels';

            for o = 1:length(sprite_types)
                for level = all_levels
                    which = strcmp(game_names, game_name) & strcmp(object_types, sprite_types{o}) & levels == level;
                    interactions{g, a}{o}(s, level) = sum(which); % {game, agent}{object type}(subject, level) = # interactions
                end
            end
        end
    end
end

%% plot by game

figure('pos', [64 221 683 1038]);

all_ms = [];

for g = 1:length(all_game_names)
    game_name = all_game_names{g};
    sprite_types = sprite_types_ordered{g};
    
    clear sems
    clear ms;
    for a = 1:length(agents)
        agent_name = agents(a).name;

        D = []; % subj x object_types
        for o = 1:length(sprite_types)
            D(:, o) = nanmean(interactions{g, a}{o}, 2); % average across levels, to get subject-level interaction counts
        end
        if ~count_walls
            D(:, strcmp(sprite_types, 'wall')) = NaN; % omit walls
        end
        [sems(a, :), ms(a, :)] = wse(D);
    end
    all_ms = [all_ms, ms];

    subplot(length(all_game_names), 1, g);

    h = bar(ms);
    hold on;

    % do some shenanigans to get the x ticks of the bars
    xs = [];
    for i = 1:length(h)
        xs = [xs, h(i).XData + h(i).XOffset];
    end
    xs = sort(xs);
    ms = ms'; ms = ms(:);
    sems = sems'; sems = sems(:); 

    errorbar(xs, ms, sems, 'o', 'MarkerSize', 1, 'color', 'black');

    % dividers
    ax = gca;
    for a = 1:length(agents)-1
        plot([a + 0.5, a + 0.5], ax.YLim, '--', 'color', [0.3 0.3 0.3]);
    end

    ylabel('# interactions / level');
    xticklabels({agents.name});
    legend(sprite_types);
    title(game_name, 'interpreter', 'none');

    hold off;
end


%% plot by level

for a = 1:length(agents)

    figure('pos', [64 421 2282 838]);

    agent_name = agents(a).name;

    for g = 1:length(all_game_names)
        game_name = all_game_names{g};
        sprite_types = sprite_types_ordered{g};
        
        clear sems
        clear ms;
        for level = all_levels
            D = []; % subj x object_types
            for o = 1:length(sprite_types)
                D(:, o) = interactions{g, a}{o}(:, level);
            end
            if ~count_walls
                D(:, strcmp(sprite_types, 'wall')) = NaN; % omit walls
            end
            [sems(level, :), ms(level, :)] = wse(D);
        end

        subplot(length(all_game_names), 1, g);

        h = bar(ms);
        hold on;

        % do some shenanigans to get the x ticks of the bars
        xs = [];
        for i = 1:length(h)
            xs = [xs, h(i).XData + h(i).XOffset];
        end
        xs = sort(xs);
        ms = ms'; ms = ms(:);
        sems = sems'; sems = sems(:); 

        errorbar(xs, ms, sems, 'o', 'MarkerSize', 1, 'color', 'black');
        ylabel('# interactions / level');
        xticks(all_levels);
        xlabel('level');
        legend(sprite_types);
        title(game_name, 'interpreter', 'none');

        hold off;
    end

end
