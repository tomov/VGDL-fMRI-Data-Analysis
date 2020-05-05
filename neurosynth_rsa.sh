mkdir output

outfileprefix="output/neurosynth_rsa"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- Running neurosynth_rsa  >> jobs.txt
echo ---------------- >> jobs.txt


# function neurosynth_rsa(rsa_idx, use_smooth, lateralized, nperms, roi_idx_min, roi_idx_max, subbatch_size)

declare -a fn_calls=(
                     "neurosynth_rsa(2, true, true, 100, 1, 400, 50)"
                     "neurosynth_rsa(1, true, true, 100, 1, 400, 50)"
                     "neurosynth_rsa(1, false, true, 100, 1, 400, 50)"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf --mem 10001 -t 4-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo neurosynth_rsa.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

