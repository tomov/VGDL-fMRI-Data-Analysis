mkdir output

outfileprefix="output/neurosynth_cva"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running neurosynth_cva  >> jobs.txt
echo ---------------- >> jobs.txt


# function neurosynth_cva(cva_idx, use_smooth, lateralized, nperms, parcel_idx, subbatch_size)

declare -a fn_calls=(
                     "neurosynth_cva(1, false, true, [], {\'vgfmri3_bait\', \'vgfmri3_chase\', \'vgfmri3_helper\', \'vgfmri3_lemmings\', \'vgfmri3_plaqueAttack\', \'vgfmri3_zelda\'}, 0.05)"
                     )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p ncf --mem 20001 -t 4-1:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo neurosynth_cva.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

