mkdir output

outfileprefix="output/might"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running might  >> jobs.txt
echo ---------------- >> jobs.txt


# function function might(rsa_idx, r, use_smooth, method, zsc)

declare -a fn_calls=(
                     "might(1, 4, false, 'lda_shrinkage', 'none')"
                     "might(1, 10, false, 'lda_shrinkage', 'none')"
                     "might(5, 4, false, 'lda_shrinkage', 'none')"
                     "might(5, 10, false, 'lda_shrinkage', 'none')"
                     "might(6, 4, false, 'lda_shrinkage', 'none')"
                     "might(6, 10, false, 'lda_shrinkage', 'none')"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf --mem 10001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo might.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

