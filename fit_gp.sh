mkdir output

#subjects=( 1 2 3 4 5 6 7 8 )  #  e.g. subjects=( 1 2 5 6 7 10 )
subjects=( 8  )  #  e.g. subjects=( 1 2 5 6 7 10 )
subj_arg="${subjects[@]}" # stringify it

what="sprite"
mask="masks/mask.nii"
glmodel=21
use_smooth=true

echo ---------------- >> jobs.txt
echo --- $(date): Running fit_gp for subjects ${subj_arg} in parallel >> jobs.txt
echo ---------------- >> jobs.txt
head -n 1 gitlog.txt >> jobs.txt

for subj in ${subjects[*]}; do
        outfileprefix="output/fit_gp_${subj}_${use_smooth}_${glmodel}_${what}"
        echo ---------------------------------------------------------------------------------
        echo Subject ${subj}, file prefix = $outfileprefix

        # send the job to NCF
        #
        sbatch_output=`sbatch -p ncf --mem 20001 -t 2-12:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'fit_gp(${subj}, ${use_smooth}, ${glmodel}, \'${mask}\', \'${what}\');exit'"`
        # for local testing
        #sbatch_output=`echo Submitted batch job 88725418`
        echo $sbatch_output

        # Append job id to jobs.txt
        #
        sbatch_output_split=($sbatch_output)
        job_id=${sbatch_output_split[3]}
        echo fit_gp.sh for subject ${subj}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

        echo watch job status with: sacct -j ${job_id}
        echo watch output with: tail -f ${outfileprefix}_${job_id}.out
        echo watch error with: tail -f ${outfileprefix}_${job_id}.err

        sleep 1
done
