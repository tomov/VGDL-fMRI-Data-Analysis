# run ccnl_fmri_preproc for a bunch of subjects
#

mkdir output

subjects=( 8  )  #  e.g. subjects=( 1 2 5 6 7 10 )
subj_arg="${subjects[@]}" # stringify it

echo ---------------- >> jobs.txt
echo --- $(date): Running ccnl_fmri_preproc for subjects ${subj_arg} in parallel >> jobs.txt
echo ---------------- >> jobs.txt

for subj in ${subjects[*]}; do
    outfileprefix="output/ccnl_fmri_preproc_${subj}"
    echo ---------------------------------------------------------------------------------
    echo Subject ${subj}, file prefix = $outfileprefix

    # send the job to NCF
    #
    sbatch_output=`sbatch -p ncf --mem 50000 -t 1-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'ccnl_fmri_preproc(vgdl_expt(),  [$subj]);exit'"`
    # for local testing
    #sbatch_output=`echo Submitted batch job 88725418`
    echo $sbatch_output

    # Append job id to jobs.txt
    #
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo ccnl_fmri_preproc.sh for subject ${subj}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    echo watch job status with: sacct -j ${job_id}
    echo watch output with: tail -f ${outfileprefix}_${job_id}.out
    echo watch error with: tail -f ${outfileprefix}_${job_id}.err

    sleep 1
done
