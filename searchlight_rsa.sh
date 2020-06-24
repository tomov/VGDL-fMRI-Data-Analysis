mkdir output

outfileprefix="output/searchlight_rsa"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running searchlight_rsa  >> jobs.txt
echo ---------------- >> jobs.txt


# function searchlight_rsa(rsa_idx, use_smooth, lateralized, nperms, parcel_idx, subbatch_size)

declare -a fn_calls=(
                     "searchlight_rsa(1, true, true, 0, 500, 6)"
                     "searchlight_rsa(1, false, true, 0, 500, 6)"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf --mem 20001 -t 1-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo searchlight_rsa.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

