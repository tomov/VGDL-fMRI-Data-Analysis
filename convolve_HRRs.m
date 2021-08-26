%
% from convolve_HRRs() in HRR.py
%
function [Xx, r_id] = convolve_HRRs(HRRs, ts, run_id, SPM)


    %nruns = length(SPM.nscan);
    nruns = double(max(run_id));
    assert(nruns <= 6);

    Xx = [];
    r_id = [];

    for s = 1:nruns
        
        which = run_id == s;

        % from spm_get_ons.m
        %
        k = SPM.nscan(s);
        assert(k == 283);

        T = SPM.xBF.T;
        assert(T == 16);

        dt = SPM.xBF.dt;
        assert(dt == 0.1250);

        UNITS = SPM.xBF.UNITS;
        assert(isequal(SPM.xBF.UNITS, 'secs'));
        TR = 1;

        bf = SPM.xBF.bf;
        assert(size(bf,1) == 257);

        % from spm_get_ons.m
        %

        ons = ts(which);
        u = HRRs(which,:);
        ton = round(ons*TR/dt) + 33; % 32 bin offset
        %toff = ton + 1; % for sanity, with impulse regressors only, e.g. theory_change_flag w/ GLM 3
        toff = ton(2:end); % frames are back-to-back
        toff = [toff; ton(end) + T]; % 1 s duration for last one, to match other between-block / level durations TODO more rigorous
        sf = zeros((k*T + 128), size(u,2));

        assert(all(ton >= 0));
        assert(all(ton < size(sf,1)));
        assert(all(toff >= 0));
        assert(all(toff < size(sf,1)));

        for j = 1:length(ton)
            sf(ton(j),:) = sf(ton(j),:) + u(j,:);
            sf(toff(j),:) = sf(toff(j),:) - u(j,:);
        end
        
        sf = cumsum(sf);
        sf = sf(1:(k*T + 32),:);                 %  32 bin offset

        % from spm_Volterra.m
        %

        %{
        % don't use convolution matrix; it's slower...
        if ~exist('A', 'var')
            A = convmtx(bf, size(sf,1));
            d = 1:size(sf,1);
        end
        X = A * sf;
        X = X(d,:);
        %}
        X = zeros(size(sf));
        for i = 1:size(sf,2)
            x = sf(:,i);
            d = 1:length(x);
            x = conv(x, bf);
            x = x(d);
            X(:,i) = x;
        end

        % from spm_fMRI_design.m
        %

        fMRI_T = SPM.xBF.T;
        fMRI_T0 = SPM.xBF.T0;

        assert(k == 283);
        idx = (0:(k - 1))*fMRI_T + fMRI_T0 + 32;
        X = X(idx,:);

        Xx = [Xx; X];
        r_id = [r_id; ones(size(X,1),1) * s];
    end

end


