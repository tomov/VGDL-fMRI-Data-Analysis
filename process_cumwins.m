function data = process_cumwins(cumwins, max_human_steps)
    % in: [games x subjects x time steps]
    % out: [subjects x time steps], time steps cropped off
    %
    data = mean(cumwins,1);
    %data = squeeze(data); % doesn't work if there is only one subject
    data = permute(data, [2,3,1]);
    data = data(:,1:max_human_steps);
end
