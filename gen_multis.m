% generate multi.mat files, for NCF (which doesn't support Mongo)

glmodels = [129:132]

EXPT = vgdl_expt;

for glmodel = glmodels
    for subj = 1:length(EXPT.subject)
        for run = 1:length(EXPT.subject(subj).functional)
            %multi = vgdl_create_multi(glmodel, subj, run)
            ccnl_check_multi(EXPT, glmodel, subj, run)
        end
    end
end
