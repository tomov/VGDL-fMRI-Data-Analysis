mkdir output

outfileprefix="output/confirmatory_betas_for_masks"
echo File prefix = $outfileprefix

echo ---------------- >> jobs.txt
echo --- $(date) Running confirmatory_betas_for_masks  >> jobs.txt
echo ---------------- >> jobs.txt


# function function confirmatory_betas_for_masks(rsa_idx, r, use_smooth, method, zsc)

declare -a fn_calls=(
                     "confirmatory_betas_for_masks(21, \'theory_change_flag\', 1, 4)"
                     "confirmatory_betas_for_masks(21, \'theory_change_flag\', 1, 6)"
                     "confirmatory_betas_for_masks(21, \'theory_change_flag\', 1, 10)"
                     "confirmatory_betas_for_masks(21, \'theory_change_flag\', 3, 4)"
                     "confirmatory_betas_for_masks(21, \'theory_change_flag\', 3, 6)"
                     "confirmatory_betas_for_masks(21, \'theory_change_flag\', 3, 10)"
                     )

declare -a fn_calls=(
                     "confirmatory_betas_for_masks(\'AAL2_GP_EMPA_grouped\', \'\', 0, 0)"
                     "confirmatory_betas_for_masks(\'AAL2_GP_EMPA\', \'\', 0, 0)"
                    )

declare -a fn_calls=(
                     "confirmatory_betas_for_masks(\'AAL2_GLM_102_grouped\', \'\', 0, 0)"
                     "confirmatory_betas_for_masks(\'AAL2_GLM_102\', \'\', 0, 0)"
                    )


for fn_call in "${fn_calls[@]}"
do
    echo $fn_call
    sbatch_output=`sbatch -p fasse --mem 10001 -t 0-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'${fn_call};exit'"`
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo $sbatch_output
    echo confirmatory_betas_for_masks.sh $fn_call: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    sleep 1
done

