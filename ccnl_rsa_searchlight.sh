# run ccnl_rsa_searchlight in batches for a range of rsa_idx
#

mkdir output

goodSubjects=( 1 2 3 4 5 6 7 8 )  # same as getGoodSubjects(), e.g. goodSubjects = ( 1 2 3 5 7 10 )
#goodSubjects=( 1 2  )  # same as getGoodSubjects(), e.g. goodSubjects = ( 1 2 3 5 7 10 )
subj_arg="${goodSubjects[@]}" # stringify it

echo ---------------- >> jobs.txt
echo --- $(date): Running ccnl_rsa_searchlight for subjects ${subj_arg} >> jobs.txt
echo ---------------- >> jobs.txt

batch_size=10000;
subbatch_size=1000;

for rsa_idx in {1..1}
do
    echo -- rsa_idx ${rsa_idx} -- >> jobs.txt

    for batch in {1..2} # 1..25
    do
        start_idx=$(((batch - 1) * batch_size + 1))
        end_idx=$((batch * batch_size))

        shuffledSubjects=( $(printf '%s\n' "${goodSubjects[@]}" | shuf ) )   # shuffle subjects so parallel GLM's don't use the same hard disk
        subj_arg="${shuffledSubjects[@]}" # stringify it

        outfileprefix="output/ccnl_rsa_searchlight_${rsa_idx}_${start_idx}-${end_idx}_goodSubjects"
        echo File prefix = $outfileprefix

        rsa_searchlight_call="ccnl_rsa_searchlight(vgdl_expt(), $rsa_idx, ${start_idx}:${end_idx}, $subbatch_size, [$subj_arg])"
        echo $rsa_searchlight_call

        # send the job to NCF
        #
        sbatch_output=`sbatch -p ncf --mem 50000 -t 4-18:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $\"$rsa_searchlight_call;exit\""`
        # for local testing
        #sbatch_output=`echo Submitted batch job 88725418`
        echo $sbatch_output

        # Append job id to jobs.txt
        #
        sbatch_output_split=($sbatch_output)
        job_id=${sbatch_output_split[3]}
        echo ccnl_rsa_searchlight.sh for rsa_idx ${rsa_idx}, ${start_idx}-${end_idx}, subjects ${subj_arg}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

        echo watch job status with: sacct -j ${job_id}
        echo watch output with: tail -f ${outfileprefix}_${job_id}.out
        echo watch error with: tail -f ${outfileprefix}_${job_id}.err

        sleep 1
    done
done
