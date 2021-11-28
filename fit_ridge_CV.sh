mkdir output

#subjects=( 1 2 3 4 5 6 7 8 )  #  e.g. subjects=( 1 2 5 6 7 10 )
#subjects=( 1  )  #  e.g. subjects=( 1 2 5 6 7 10 )
subjects=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31  32  )  #  e.g. subjects=( 1 2 5 6 7 10 )
#subjects=( 2 3 4 5 6 7 8 )  #  e.g. subjects=( 1 2 5 6 7 10 )
subj_arg="${subjects[@]}" # stringify it

#model_name="DQN"
#model_name="EMPA"
model_name="state"
#what="conv1"
#what="theory"
what=""
mask="masks/mask.nii"
glmodel=1
use_smooth=true
subsample_only=false
project=true
save_Y_hat=false

echo ---------------- >> jobs.txt
echo --- $(date): Running fit_ridge_CV for subjects ${subj_arg} in parallel >> jobs.txt
echo ---------------- >> jobs.txt
head -n 1 gitlog.txt >> jobs.txt

for subj in ${subjects[*]}; do
    outfileprefix="output/fit_ridge_CV_${subj}_${use_smooth}_${glmodel}"
    echo ---------------------------------------------------------------------------------
    echo Subject ${subj}, file prefix = $outfileprefix

    # send the job to NCF
    #
    sbatch_output=`sbatch -p  fasse --mem 20001 -t 0-15:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'fit_ridge_CV(${subj}, ${use_smooth}, ${glmodel}, \'${mask}\', \'${model_name}\', \'${what}\', ${subsample_only}, ${project}, ${save_Y_hat});exit'"`
    # for local testing
    #sbatch_output=`echo Submitted batch job 88725418`
    echo $sbatch_output

    # Append job id to jobs.txt
    #
    sbatch_output_split=($sbatch_output)
    job_id=${sbatch_output_split[3]}
    echo fit_ridge_CV.sh for subject ${subj}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

    echo watch job status with: sacct -j ${job_id}
    echo watch output with: tail -f ${outfileprefix}_${job_id}.out
    echo watch error with: tail -f ${outfileprefix}_${job_id}.err

    sleep 1
done
