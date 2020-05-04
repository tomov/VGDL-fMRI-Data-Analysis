function EXPT = vgdl_expt_nosmooth(local)

    % same as vgdl_expt but for nonsmoothed data
    % call with ccnl_fmri_glm_nosmooth

    EXPT = vgdl_expt();

    EXPT.modeldir = [EXPT.exptdir, 'glmOutput_nosmooth'];
    EXPT.rsadir = [EXPT.rsadir, 'rsaOutput_nosmooth'];

