# run ccnl_fmri_con for a bunch of subjects; must edit run_ccnl_fmri_con.m first (that's where the action is)
#

mkdir output

goodSubjects=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 )  # same as getGoodSubjects(), e.g. goodSubjects = ( 1 2 3 5 7 10 )
subj_arg="${goodSubjects[@]}" # stringify it

models=(55 56 51 52 71 79 85 )

echo ---------------- >> jobs.txt
echo --- $(date): Running ccnl_fmri_con for GLMs ${models}, subjects ${subj_arg} >> jobs.txt
echo ---------------- >> jobs.txt


for model in "${models[@]}"
do
    shuffledSubjects=( $(printf '%s\n' "${goodSubjects[@]}" | shuf ) )   # shuffle subjects so parallel GLM's don't use the same hard disk
    subj_arg="${shuffledSubjects[@]}" # stringify it

    outfileprefix="output/ccnl_fmri_con_${model}_goodSubjects"
    echo File prefix = $outfileprefix

    # send the job to NCF
    #
    sbatch_output=`sbatch -p ncf --mem 5001 -t 0-12:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'run_ccnl_fmri_con($model, [$subj_arg]);exit'"`
    # for local testing
    #sbatch_output=`echo Submitted batch job 88725418`
    echo $sbatch_output

    # Append job id to jobs.txt
    #
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo ccnl_fmri_con.sh for GLM ${model}, subjects ${subj_arg}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    echo watch job status with: sacct -j ${job_id}
    echo watch output with: tail -f ${outfileprefix}_${job_id}.out
    echo watch error with: tail -f ${outfileprefix}_${job_id}.err

    sleep 1
done

