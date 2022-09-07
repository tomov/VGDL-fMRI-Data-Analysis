# run ccnl_fmri_glm (single subject GLM) for a range of models for a bunch of subjects
# copied from Exploration project https://github.com/tomov/Exploration-Data-Analysis
#

mkdir output

# this is crucial -- we can't simulate the subjects that don't have the full data; we don't want ccnl_fmri_glm to error out (note that when we were doing them in parallel before, it didn't matter if one subject failed b/c all other jobs would still continue)
#
goodSubjects=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 )  # same as getGoodSubjects(), e.g. goodSubjects = ( 1 2 3 5 7 10 )
subj_arg="${goodSubjects[@]}" # stringify it

#models=( 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 )
#models=( 166  168 170 172 174 176 178 180 )
models=( 196 )
#models=( 181 182 183 184 )

echo ---------------- >> jobs.txt
echo --- $(date): Running ccnl_fmri_glm for subjects ${subj_arg} >> jobs.txt
echo ---------------- >> jobs.txt


for model in "${models[@]}"
do
    shuffledSubjects=( $(printf '%s\n' "${goodSubjects[@]}" | shuf ) )   # shuffle subjects so parallel GLM's don't use the same hard disk
    subj_arg="${shuffledSubjects[@]}" # stringify it

    outfileprefix="output/ccnl_fmri_glm_${model}_goodSubjects"
    echo File prefix = $outfileprefix

    # send the job to NCF
    #
    sbatch_output=`sbatch -p fasse --mem 10000 -t 1-12:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'ccnl_fmri_glm(vgdl_expt(), $model, [$subj_arg]);exit'"`
    # for local testing
    #sbatch_output=`echo Submitted batch job 88725418`
    echo $sbatch_output

    # Append job id to jobs.txt
    #
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo ccnl_fmri_glm.sh for GLM ${model}, subjects ${subj_arg}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    echo watch job status with: sacct -j ${job_id}
    echo watch output with: tail -f ${outfileprefix}_${job_id}.out
    echo watch error with: tail -f ${outfileprefix}_${job_id}.err

    sleep 1
done
