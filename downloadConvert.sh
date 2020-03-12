#!/bin/bash
#SBATCH -p ncf # partition (queue)
#SBATCH --mem 8000 # memory
#SBATCH -t 0-8:00 # time (D-HH:MM)

# paremeter = subject ID on CBS central, e.g. 200311_VGDL_001
# USAGE ex: ./downloadConvert.sh 200311_VGDL_001

# WHEN USING ON NEW STUDY:
#
# 1) look for 'experiment' and make sure dinectory structure for data and scripts conforms the format
# 2) make sure data is on CBS central
# 3) make sure ArcGet.py script is in same dir
# 4) make sure subj id is correct

# download all structural and functional scans in order in which they were taken
#
python ArcGet.py -a cbscentral -s ${1} -r MEMPRAGE\ RMS,Minn_HCP_1.7mm_S3p2_Task

experiment="VGDL_fMRI" # must match directory name

# fileNames should correspond to order in which scans were taken, e.g. normally it's fileNames=(struct run001 run002 run003 run004 run005 run006 run007 run008)
# HOWEVER, if things are out of order or some scans were bad, might have to reorder and include dummy entries for the bad scans,
# e.g. if the structural was bad and you took another structural in the end, it would look like this: fileNames=(struct_bad run001 run002 run003 run004 run005 run006 run007 run008 struct)
fileNames=(struct run001 run002 run003 run004 run005 run006)


mkdir /ncf/gershman/Lab/${experiment}/subjects/
cd /ncf/gershman/Lab/${experiment}/subjects/

mkdir /ncf/gershman/Lab/${experiment}/subjects/${1}/
mkdir  /ncf/gershman/Lab/${experiment}/subjects/${1}/RAW
mkdir /ncf/gershman/Lab/${experiment}/subjects/${1}/preproc

#ArcGet.py -a cbscentral -s ${1} 

cd /ncf/gershman/Lab/${experiment}/subjects/${1}/RAW

myruns=`ls *.MR.Investigators_Gershman.*.1.* | sort -t. -nk4`

count=0
for run in ${myruns[*]}; do
	echo ------------------
	echo Run = ${run}
	echo count = ${count}
	echo filename = ${fileNames[$count]}
	mri_convert -it siemens_dicom -ot nii -i ${run} -o /ncf/gershman/Lab/${experiment}/subjects/${1}/preproc/${fileNames[$count]}.nii
	count=${count}+1
done
