function substring = get_unique_regressor_substring(name)
    if contains(name, ') vgfmri')
        % lump together all games under the same regressor
        substring = ') vgfmri';
    else
        % just get the regressor name without the session
        substring = name(5:end);
    end
