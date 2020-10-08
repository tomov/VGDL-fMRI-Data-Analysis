function loglik = lik_empa(x, data)

    % Eps-greedy likelihood f'n for EMPA based on avatar-sprite interactions
    %
    % INPUTS:
    %   x - parameters:
    %     x(1) - eps for eps-greedy 
    %   data - struct:
    %     behavior - {N x 1} cell array of sprites avatar touched (i.e. interacted with)
    %     predictions - {N x 1} cell array of sprites EMPA predicted avatar will touch next after previous interaction
    %
    % OUTPUTS:
    %   loglik - log-likelihood

    eps = x(1);

    B = data.behavior;
    P = data.predictions;

    assert(length(B) == length(P));

    loglik = 0;
    for n = 1:length(B)
        b = B{n}; % subject interactions
        b = arrayfun(@strtrim, mat2cell(b, ones(1,size(b,1)))); % format into cell array
        b(contains(b, 'DARKGRAY')) = []; % remove walls
        if isempty(b)
            continue
        end
        
        p = P{n};
        p = arrayfun(@strtrim, mat2cell(p, ones(1,size(p,1))));
        if isempty(p)
            continue % TODO should not happen in theory
        end
        p(contains(p, 'DARKGRAY')) = []; % remove walls

        is_pred = any(ismember(b, p)); % TODO too lenient
        if is_pred
            % subject interacted with EMPA-predicted sprite
            loglik = loglik + log(1 - eps);
        else 
            loglik = loglik + log(eps);
        end
    end

