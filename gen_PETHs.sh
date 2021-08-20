mkdir output

outfileprefix="output/gen_PETHs"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running gen_PETHs  >> jobs.txt
echo ---------------- >> jobs.txt


# function function gen_PETHs(rsa_idx, r, use_smooth, method, zsc)

declare -a fn_calls=(
                     "gen_PETHs(21, \'theory_change_flag\', 1, 4)"
                     "gen_PETHs(21, \'theory_change_flag\', 3, 4)"
                     "gen_PETHs(21, \'theory_change_flag\', 3, 10)"
                     "gen_PETHs(\'tomov2018KL\', [], [], 4)"
                     "gen_PETHs(\'tomov2018KL\', [], [], 10)"
                     "gen_PETHs(\'hayley2021psi\', [], [], 4)"
                     "gen_PETHs(\'hayley2021psi\', [], [], 10)"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf --mem 10001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo gen_PETHs.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

