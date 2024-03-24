# purpose - merge biological replicates and already paired technical replicates (now treated as biological replicates) of ChIP libraries
# usage
# qq
#
# example
#
# 'qq'

#input directory
INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/all_reps_dedup

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=10:00:00
#PBS -A open
cd $OUTPUT

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

	#Get basename for PBS
	#set output file names
	FINAL_MERGE=$(echo "K562_benzonase-seq_final_merge.bam")
	MASTER=$(echo "K562_benzonase-seq_master.bam")

	sampleID=master\_merge_filter.pbs
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output directory" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#merge all deduplicated replicates" >> $sampleID
	echo "samtools merge $FINAL_MERGE $INPUT/*DEDUP.bam" >> $sampleID
	echo "#filter bam of all merged replicates" >> $sampleID
	echo "samtools view -o $MASTER -h -b -f 0x1 -F 0x404 $FINAL_MERGE" >> $sampleID
	echo "#remove intermediate file" >> $sampleID
	echo "rm $FINAL_MERGE" >> $sampleID
	echo "#finish script" >> $sampleID
