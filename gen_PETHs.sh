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

declare -a fn_calls=(
                     "gen_PETHs(\'AAL2_GP_EMPA_grouped\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_grouped\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GP_EMPA\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL2_GP_EMPA\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GLM_102_grouped\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL2_GLM_102_grouped\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GLM_102\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL2_GLM_102\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP\')"
                     )

# Figure 5
declare -a fn_calls=(
                     "gen_PETHs(\'AAL2_GP_EMPA_grouped\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GP_EMPA\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GLM_102_grouped\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GLM_102\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP\')"
                     )

# Figure 5 ++
#declare -a fn_calls=(
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_sprite\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_sprite\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_interaction\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_interaction\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_termination\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_termination\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_sprite\', 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_sprite\', 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_interaction\', 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_interaction\', 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_termination\', 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_termination\', 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_sprite\', 1, 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_sprite\', 1, 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_interaction\', 1, 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_interaction\', 1, 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_termination\', 1, 1)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_termination\', 1, 1)"
#                     )

#declare -a fn_calls=(
#                     "gen_PETHs(\'AAL2_GP_EMPA_grouped\', \'\', 0, 0, \'GP_DQN\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA\', \'\', 0, 0, \'GP_DQN\')"
#                     "gen_PETHs(\'AAL2_GLM_102_grouped\', \'\', 0, 0, \'GP_DQN\')"
#                     "gen_PETHs(\'AAL2_GLM_102\', \'\', 0, 0, \'GP_DQN\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_DQN\')"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_DQN\')"
#                     )

declare -a fn_calls=(
                     "gen_PETHs(\'AAL2_GP_EMPA_grouped\', \'\', 0, 0, \'GP_VAE\')"
                     "gen_PETHs(\'AAL2_GP_EMPA\', \'\', 0, 0, \'GP_VAE\')"
                     "gen_PETHs(\'AAL2_GLM_102_grouped\', \'\', 0, 0, \'GP_VAE\')"
                     "gen_PETHs(\'AAL2_GLM_102\', \'\', 0, 0, \'GP_VAE\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'GP_VAE\')"
                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'GP_VAE\')"
                     )

# Figure 4 without baseline
#
#declare -a fn_calls=(
#                     "gen_PETHs(\'AAL2_GP_EMPA_grouped\', \'\', 0, 0, \'BOLD\', false, true)"
#                     "gen_PETHs(\'AAL2_GP_EMPA\', \'\', 0, 0, \'BOLD\\', false, true)"
#                     "gen_PETHs(\'AAL2_GLM_102_grouped\', \'\', 0, 0, \'BOLD\\', false, true)"
#                     "gen_PETHs(\'AAL2_GLM_102\', \'\', 0, 0, \'BOLD\\', false, true)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102_grouped\', \'\', 0, 0, \'BOLD\\', false, true)"
#                     "gen_PETHs(\'AAL2_GP_EMPA_GLM_102\', \'\', 0, 0, \'BOLD\\', false, true)"
#                     )

# your on
# Neuron R1

declare -a fn_calls=(
                     "gen_PETHs(\'Brodmann\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'Brodmann\', \'\', 0, 0, \'GP\')"
                     "gen_PETHs(\'AAL3v1_neuron\', \'\', 0, 0, \'BOLD\')"
                     "gen_PETHs(\'AAL3v1_neuron\', \'\', 0, 0, \'GP\')"
                     )

declare -a fn_calls=(
                     "gen_PETHs(\'Brodmann\', \'\', 0, 0, \'GP_sprite\')"
                     "gen_PETHs(\'Brodmann\', \'\', 0, 0, \'GP_interaction\')"
                     "gen_PETHs(\'Brodmann\', \'\', 0, 0, \'GP_termination\')"
                     "gen_PETHs(\'AAL3v1_neuron\', \'\', 0, 0, \'GP_sprite\')"
                     "gen_PETHs(\'AAL3v1_neuron\', \'\', 0, 0, \'GP_interaction\')"
                     "gen_PETHs(\'AAL3v1_neuron\', \'\', 0, 0, \'GP_termination\')"
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

