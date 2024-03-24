# purpose - merge technical teplicates of benzonase-seq
# usage
# qq
#
# example
#
# 'qq'

#technical replicates
REP1=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/raw_aligned/33028_Input_-_K562_raw.bam
REP2=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/raw_aligned/33164_Input_-_K562_raw.bam

#input directory
INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/raw_aligned

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/tech_reps_merged

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=3:00:00
#PBS -A open
cd $OUTPUT

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

for file in $INPUT; do
	#Get basename for PBS
	fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	#set output file names
	REP1a=$(echo $REP1 | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | awk '{print $1}')
	REP2a=$(echo $REP2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	MERGE=$(echo "$REP1a""_""$REP2a" | awk -F. '{print $1"_MERGE.bam"}')
	MERGEa=$(echo "$REP1a""_""$REP2a" | awk -F. '{print $1"_MERGE"}')

	sampleID=$MERGEa\.pbs
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output directory" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#merge technical replicates. input files are already sorted" >> $sampleID
	echo "samtools merge $MERGE $REP1 $REP2" >> $sampleID
	echo "# finish script" >> $sampleID

done
