%
% like convolve_HRRs but subsampling only, no convolution
%
function [Xx, r_id] = subsample_HRRs(HRRs, ts, run_id, SPM)

    TR = 2; % TODO hardcoded

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

        % from spm_get_ons.m
        %

        ons = ts(which,:);
        u = HRRs(which,:);
       
        t = 1:TR:k*TR;
        [~,ix] = min(abs(ons - t), [], 1);

        X = u(ix,:);
        Xx = [Xx; X];
        r_id = [r_id; ones(size(X,1),1) * s];
    end

end


