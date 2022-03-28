mkdir output

#subjects=( 1 2 3 4 5 6 7 8 )  #  e.g. subjects=( 1 2 5 6 7 10 )
#subjects=( 1  )  #  e.g. subjects=( 1 2 5 6 7 10 )
subjects=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  32  )  #  e.g. subjects=( 1 2 5 6 7 10 )
#subjects=( 16 )  #  e.g. subjects=( 1 2 5 6 7 10 )
#subjects=(  17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  32  )
subj_arg="${subjects[@]}" # stringify it

mask="masks/mask.nii"
#model_name="EMPA"
#model_name="game"
#model_name="irrelevant"
#model_name="state"
#model_name="DQN"
model_name="DQN25M"
#model_name="PCA"
#model_name="VAE"
#what="conv3"
what="conv1"
#what="all"
#what="theory"
#what=""
#what="sprite"
#what="interaction"
#what="termination"
#what="novelty"
glmodel=1
use_smooth=true
project=1
normalize=1
concat=0
novelty=1
fast=true
save_Y_hat=0

echo ---------------- >> jobs.txt
echo --- $(date): Running fit_gp_CV for subjects ${subj_arg} in parallel >> jobs.txt
echo ---------------- >> jobs.txt
head -n 1 gitlog.txt >> jobs.txt

for subj in ${subjects[*]}; do
    outfileprefix="output/fit_gp_CV_${subj}_${use_smooth}_${glmodel}_${model_name}_${what}_${project}_${normalize}_${save_Y_hat}"
    echo ---------------------------------------------------------------------------------
    echo Subject ${subj}, file prefix = $outfileprefix

    # send the job to NCF
    #
    sbatch_output=`sbatch -p fasse --mem 20001 -t 0-5:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'fit_gp_CV(${subj}, ${use_smooth}, ${glmodel}, \'${mask}\', \'${model_name}\', \'${what}\', ${project}, ${normalize}, ${concat}, ${novelty}, ${fast}, ${save_Y_hat});exit'"`
    # for local testing
    #sbatch_output=`echo Submitted batch job 88725418`
    echo $sbatch_output

    # Append job id to jobs.txt
    #
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo fit_gp_CV.sh for subject ${subj}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    echo watch job status with: sacct -j ${job_id}
    echo watch output with: tail -f ${outfileprefix}_${job_id}.out
    echo watch error with: tail -f ${outfileprefix}_${job_id}.err

    sleep 1
done
