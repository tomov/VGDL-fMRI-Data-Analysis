% generate rsa.mat files, for NCF (which doesn't support Mongo)
% copy of gen_multis.m

rsa_idxs = [3]

EXPT = vgdl_expt;

for glmodel = glmodels
    for subj = 1:length(EXPT.subject)
        for run = 1:length(EXPT.subject(subj).functional)
            rsa = vgdl_create_rsa(rsa_idxs, subj)
        end
    end
end
