# download fMRI data for a bunch of subjects using ArcGet and convert it to nice .nii format using mri_convert
#

module load mri_convert/2015_12_03-ncf

mkdir output

experiment="VGDL_fMRI" # must match directory name TODO IMPORTANT script directory hardcoded below

# list of subject ID's as registered on CBS, e.g.  subjects=('180725_UEP_001' '189725_UEP_002')
subjects=( '200315_VGDL_008')
#subjects=('180725_UEP_001' '189725_UEP_002')

#subjects=$(/ncf/gershman/Lab/${experiment}/subjects.txt)


echo ---------------- >> jobs.txt
echo --- Running downloadConvertSBatch for subjects ${subjects} >> /ncf/gershman/Lab/scripts/matlab/${experiment}/jobs.txt
echo ---------------- >> jobs.txt

for subj in ${subjects[*]}; do
    outfileprefix="/ncf/gershman/Lab/scripts/matlab/${experiment}/output/downloadConvert_${subj}"
    echo ------------------------------------------------------------------
    echo Subject ${subj}, file prefix = $outfileprefix

    # send the job to NCF
    #
    #cd /ncf/gershman/Lab/${experiment}/subjects/
    #mkdir /ncf/gershman/Lab/${experiment}/subjects/${subj}/
    sbatch_output=`sbatch -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err /ncf/gershman/Lab/scripts/matlab/${experiment}/downloadConvert.sh ${subj}` 
    # for local testing
    #sbatch_output=`echo Submitted batch job 88725418`
    echo $sbatch_output

    # Append job id to jobs.txt
    #
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo downloadConvertSBatch.sh for subject ${subj}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> /ncf/gershman/Lab/scripts/matlab/${experiment}/jobs.txt

    echo watch job status with: sacct -j ${job_id}
    echo watch output with: tail -f ${outfileprefix}_${job_id}.out
    echo watch error with: tail -f ${outfileprefix}_${job_id}.err

    sleep 1
done

