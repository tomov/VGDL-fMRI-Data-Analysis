mkdir output

outfileprefix="output/gen_PETHs"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running gen_PETHs  >> jobs.txt
echo ---------------- >> jobs.txt


# function function gen_PETHs(rsa_idx, r, use_smooth, method, zsc)

declare -a fn_calls=(
                     "gen_PETHs(\'AAL2\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL2\', \'\', 0, 0, \'GP\')"
                     )

declare -a fn_calls=(
                     "gen_PETHs(\'tomov2018KL\', 0, 0, 10)"
                     "gen_PETHs(\'hayley2021psi\', 0, 0, 10)"
                     )

declare -a fn_calls=(
                     "gen_PETHs(\'AAL3v1\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL3v1\', \'\', 0, 0, \'GP\')"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p fasse --mem 10001 -t 0-10:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo gen_PETHs.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

