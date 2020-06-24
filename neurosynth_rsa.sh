mkdir output

outfileprefix="output/neurosynth_rsa"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running neurosynth_rsa  >> jobs.txt
echo ---------------- >> jobs.txt


# function neurosynth_rsa(rsa_idx, use_smooth, lateralized, nperms, parcel_idx, subbatch_size)

declare -a fn_calls=(
                     "neurosynth_rsa(6, true, true, 10000, [177 294 365 38 293 174 194 118], 500)"
                     "neurosynth_rsa(6, false, true, 10000, [177 38 94 29 86 333 294 365 38 293 174 194 118], 500)"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf --mem 20001 -t 4-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo neurosynth_rsa.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

