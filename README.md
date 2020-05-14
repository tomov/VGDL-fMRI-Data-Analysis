# VGDL fMRI

README copied from Exploration [repo](https://github.com/tomov/Exploration-Data-Analysis).

Analyzing fMRI data of humans playing video games.

Scripts are in `/ncf/gershman/Lab/scripts/matlab/VGDL_fMRI/`.

Data are in `/ncf/gershman/Lab/VGDL_fMRI/`.

Useful links:
- [CBS cluster FAQ](http://cbs.fas.harvard.edu/science/core-facilities/neuroimaging/information-investigators/faq) -- how to use the cluster, send jobs, ArcGet.py, slurm, sacct, etc
- [CBS central login](http://cbscentral.rc.fas.harvard.edu) -- where the fMRI data live
- [ccnl_fmri wiki](https://github.com/sjgershm/ccnl-fmri/wiki) -- how to use Sam's fMRI pipeline

## To preprocess a newly scanned subject

0. Copy everything to the cluster and login
   * See `scp_to_ncf.sh`, maybe edit it accordingly and run it `./scp_to_ncf.sh`
1. Open `downloadConvertSBatch.sh` and edit `subjects` to include the **new subject only**.
2. Open `downloadConvert.sh` and edit `fileNames` if necessary. Read comments for details
3. Run `./downloadConvertSBatch.sh`
   * Takes ~30 mins to complete.
   * Make sure there are no scary errors in the .err file or the .out file.
   * Make sure all files (struct.nii, run001.nii, etc) are written normally in the .out file
   * Go to `/ncf/gershman/Lab/Exploration/subjects/` and make sure all subject files are there
4. Open `vgdl_getSubjectsDirsAndRuns.m` and **append** new subject info to `subjects`, `subjdirs`, and `nRuns` accordingly (so they include info for **all subjects**). Read comments for details.
5. Open `ccnl_fmri_preproc.sh` and edit `subjects` to include the index of the **new subject only**. This index is the ordinal of the subject info in `vgdl_getSubjectsDirsAndRuns.m`
6. Run `./ccnl_fmri_preproc.sh`
   * Takes ~10 hours (!) to complete.
   * Make sure there are no scary errors in the .err file or the .out file.
   * Make sure all files are in the data directory (compare with previous subjects).
7. Open MATLAB with a GUI and cd into scripts directory
   * Log into one of the cluster terminals in Northwest (e.g. by the fMRI scanner),
   * or download [XQuartz](https://www.xquartz.org/) on your Mac, `ssh -X` into the cluster on one of the nvdi nodes, e.g. `ssh -X mtomov13@ncfvdi.rc.fas.harvard.edu`, and run `matlab` on the command line. You should see the MATLAB splash screen
8. Run `ccnl_plot_movement(vgdl_expt(), XX)`, where XX is the subject index (e.g. 1)
   * Make sure subject didn't move too much during runs
9. Run `ccnl_check_registration(vgdl_expt(), XX)`, where XX is the subject index
   * Click around (e.g. around ventricles), make sure average functional (top) is aligned with the structural (bottom)
10. Append the subject's behavioral data to `data.csv` and upload it to the cluster

11. Check the raw BOLD traces with, e.g. `ccnl_view_mask('/ncf/gershman/Lab/VGDL_fMRI/subjects/200311_VGDL_001/preproc/swurun001.nii')`

Tips
   * Always check that the job is running with `sacct`
   * Wait until some output is printed to the .out file, to be sure that it didn't fail at the very beginning
       * particularly important for long-running jobs

Common errors and mistakes
   * ArcGet.py doesn't exist -> happens sometimes, idk why
   * ArcGet.py can't access CBS Central -> read the CBS FAQ on ArcGet.py; you need to set it up first
   * ArcGet.py can't find subject on CBS central -> make sure that the subject ID (e.g. `180807_UEP_025`) is exactly as it is on CBS Central
       * check date
       * check spaces and underscores
       * make sure no commas in the `subjects` array in the script `downloadConvertSBatch.sh`
       * make sure subject was actually sent to CBS central
       * make sure all scans were sent to CBS central
   * ArcGet.py can't write to data directory
       * make sure you have permissions or someone else didn't download that subject already
       * make sure data directory exists
       * make sure that the directories in the scripts are correct
   * Runs (run###.nii) or structurals (struct.nii) missing or empty -> make sure the file order in `fileNames` in the script `downloadConvert.sh` is the same as the session order of the subject
       * e.g. if another structural was acquired after run 4, it should be `fileNames=(struct_1 run001 run002 run004 run004 struct_2 run005 run006 run007 run008)`
       * e.g. if run 6 was interrupted and restarted, it should be `fileNames=(struct run001 run002 run004 run004 run005 run006_bad run006 run007 run008

## To create and run a new GLM

1. Add the GLM to `vgdl_create_multi.m`
2. Test the GLM locally
    * Call `multi = vgdl_create_multi(...)` in MATLAB
    * Inspect `multi` and make sure all the names, onsets, and pmods make sense
    * Batch test for all subjects and runs using `ccnl_check_multi.m` (this also creates the cached .mat files, which are important since MATLAB on the cluster doesn't have mongo)
3. Copy `vgdl_create_multi.m` to the cluster
    * See `scp_to_ncf.sh`
4. Log into the cluster and run an interactive job, e.g. with `./srun.sh`
5. Run MATLAB in the interactive job, e.g. with `./matlab.sh`
6. Test the GLM on the cluster
    * repeat same steps as testing locally
7. Dry run the GLM in MATLAB by calling `ccnl_fmri_glm(...)` to make sure SPM actually likes it; then interrupt it with Ctrl+C once it looks like it's working
    * Optionally use the `fake_glm=true` argument
8. Exit MATLAB and edit `ccnl_fmri_glm.sh`
    * Change `goodSubjects = (...)` to include all subjects
    * Change `for model in {...}` to use your GLM only
9. Run `./ccnl_fmri_glm.sh`


## To run a group-level contrast

1. Add your contrast to `run_ccnl_fmri_con.m`
    * Copy-paste one of the calls to `ccnl_fmri_con` and modify accordingly
    * Please keep all previous contrasts, just comment them out
2. Copy `run_ccnl_fmri_con.m` to the cluster
3. Log into the cluster and run an interactive job, e.g. with `./srun.sh`
4. Run MATLAB in the interactive job, e.g. with `./matlab.sh`
5. Test the contrast on the cluster
    * Call `run_ccnl_fmri_con` in MATLAB, make sure it starts doing stuff and Ctrl+C it
3. Exit MATLAB and run `./run_ccnl_fmri_con.sh`

## To visualize the group-level contrast

0. Run `ccnl_mean.m` on the cluster, scp `glmOutput/mean.nii` from cluster 
1. Run MATLAB with a GUI and cd into scripts directory
    * look for X11 above
2. Run `ccnl_view` with the corresponding GLM and contrast
    * e.g. `ccnl_view(vgdl_expt(), 1, 'RR - SS')`

## To run RSA searchlight

Preliminaries:

1. Create `masks/` directory, copy over `mask.nii` from any contrast (e.g. `../glmOutput/model1/con1/mask.nii`; it's the group-level whole-brain mask)
2. Create some anatomical masks using `ccnl_create_mask.m` -- see e.g. `create_anatomical_masks_AAL2.m
3. Create a GLM that extracts beta series -- see e.g. GLM 1 and GLM 22 for a block and event-based ways to do it
4. Run the GLM, with and without smoothing (see `vgdl_expt_nosmooth.m` and `ccnl_fmri_glm_nosmooth.sh`)

To add and run a new RSA searchlight:

5. Add a new RSA to `vgdl_create_rsa.m`
    * make sure to run e.g. `for s = 1:8, vgdl_create_rsa(2,s); end` in MATLAB to generate the `.mat` files, since MATLAB on the cluster doesn't have mongo
6. Edit `ccnl_rsa_searchlight.sh`, check and adjust accordingly: 
    * `vgdl_expt` vs. `vgdl_expt_nosmooth`,
    * subjects, 
    * batch sizes, 
    * for loop over rsa's
    * for loop over batches, 
    * memory limits,
    * time limits
7. `scp_to_ncf.sh` to copy .m and .mat files to cluster (check it first)
8. Run interactive MATLAB on cluster, test new RSA with `ccnl_check_rsa.m`
9. Run `ccnl_rsa_searchlight.sh` (recommended with only 1-2 batches first, as a sanity check)
10. `scp_from_ncf.sh` to copy `rsaOutput` files from cluster locally
11. `ccnl_rsa_view_searchlight.m` to view RSA

Note: If you re-run the beta series GLM, you will have to go inside all `rsaOutput/rsa*/` directories and remove the `beta*.mat` files; we extract and cache them to speed things up, so if the GLM changes, we have to remove them.

Note: you have to remove the .nii tmap file once it's generated by `ccnl_rsa_view` if you want to re-generate it

## To adapt pipeline to a new experiment

1. create data directory on cluster
2. create scripts directory on cluster
3. copy over `vgdl_*.m` files, rename to `yourproject_*.m`
4. create git repo, add files, push to github
5. `git grep -i 'vgdl' *.m` -- change all prefixes, directories, etc. accordingly to reflect your project

## Tips & Gotchas

Careful scp'ing to cluster while jobs are running -- reading files (e.g. mat files) can break if they're being overwritten
