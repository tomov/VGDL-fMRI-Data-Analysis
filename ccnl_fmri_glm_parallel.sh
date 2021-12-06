# run ccnl_fmri_glm (single subject GLM) for a range of models for a bunch of subjects
# copied from ../Hayley/run_ccnl_fmri_glm_parallel_jobs.sh
#

mkdir output

goodSubjects=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 )  # same as getGoodSubjects(), e.g. goodSubjects = ( 1 2 3 5 7 10 )

models=( 153 )


subj_arg="${goodSubjects[@]}" # stringify it

echo ---------------- >> jobs.txt
echo --- Running ccnl_fmri_glm_parallel for subjects ${subj_arg} >> jobs.txt
echo ---------------- >> jobs.txt


for model in "${models[@]}"
do
    shuffledSubjects=( $(printf '%s\n' "${goodSubjects[@]}" | shuf ) )   # shuffle subjects so parallel GLM's don't use the same hard disk
    subj_arg="${shuffledSubjects[@]}" # stringify it

    for subj in ${shuffledSubjects[*]}; do
        outfileprefix="output/ccnl_fmri_glm_parallel_${model}_subj_${subj}"
        echo File prefix = $outfileprefix

        # send the job to NCF
        #
        sbatch_output=`sbatch -p fasse --mem 10000 -t 0-10:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'ccnl_fmri_glm(vgdl_expt(), $model, [$subj]);exit'"`
        # for local testing
        #sbatch_output=`echo Submitted batch job 88725418`
        echo $sbatch_output

        # Append job id to jobs.txt
        #
        sbatch_output_split=($sbatch_output)
        job_id=${sbatch_output_split[3]}
        echo ccnl_fmri_glm_parallel.sh for GLM ${model}, subject ${subj}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

        echo watch job status with: sacct -j ${job_id}
        echo watch output with: tail -f ${outfileprefix}_${job_id}.out
        echo watch error with: tail -f ${outfileprefix}_${job_id}.err

        sleep 1
    done
done
