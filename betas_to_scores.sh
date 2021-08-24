mkdir output

outfileprefix="output/betas_to_scores"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running betas_to_scores  >> jobs.txt
echo ---------------- >> jobs.txt


# function function betas_to_scores(rsa_idx, r, use_smooth, method, zsc)

declare -a fn_calls=(
                     "betas_to_scores(21, \'theory_change_flag\', 1, 10, \'theory_change_flag\')"
                     "betas_to_scores(21, \'theory_change_flag\', 1, 10, \'theory_change_flag\')"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf --mem 10001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo betas_to_scores.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

