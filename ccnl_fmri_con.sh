# run ccnl_fmri_con for a bunch of subjects; must edit run_ccnl_fmri_con.m first (that's where the action is)
#

mkdir output

#goodSubjects=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 )  # same as getGoodSubjects(), e.g. goodSubjects = ( 1 2 3 5 7 10 )
goodSubjects=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14                16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 )  # EXCLUDED 15!!!!!!!!!!!!!!!!!
subj_arg="${goodSubjects[@]}" # stringify it

#models=( 27 28 29 76 82 83 81 80 86 87 88 68 30 34 31 35 36 37 42 43 46 47 48 49 57 )
#models=( 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 )
models=( 196 )
models_arg="${models[@]}" # stringify it with spaces
models_str=`echo $models_arg | sed 's/ /_/g'` # stringify with underscores

echo ---------------- >> jobs.txt
echo --- $(date): Running ccnl_fmri_con for subjects ${subj_arg} >> jobs.txt
echo ---------------- >> jobs.txt


shuffledSubjects=( $(printf '%s\n' "${goodSubjects[@]}" | shuf ) )   # shuffle subjects so parallel GLM's don't use the same hard disk
subj_arg="${shuffledSubjects[@]}" # stringify it

outfileprefix="output/ccnl_fmri_con_${models_str}_goodSubjects"
echo File prefix = $outfileprefix

# send the job to NCF
#
sbatch_output=`sbatch -p fasse --mem 10001 -t 1-12:20 -o ${outfileprefix}_%j.out -e ${outfileprefix}_%j.err --wrap="matlab -nodisplay -nosplash -nojvm -r $'run_ccnl_fmri_con([$models_arg], [$subj_arg]);exit'"`
# for local testing
#sbatch_output=`echo Submitted batch job 88725418`
echo $sbatch_output

# Append job id to jobs.txt
#
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo ccnl_fmri_con.sh for GLMs ${models_arg}, subjects ${subj_arg}: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

echo watch job status with: sacct -j ${job_id}
echo watch output with: tail -f ${outfileprefix}_${job_id}.out
echo watch error with: tail -f ${outfileprefix}_${job_id}.err
