# purpose - merge technical teplicates of benzonase-seq
# usage
# qq
#
# example
#
# 'qq'

#input directory
INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/raw_aligned

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/all_reps_dedup

#picard.jar
PICARD=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/picard.jar

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=6:00:00
#PBS -A open
cd $OUTPUT

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

for file in $INPUT/*.bam; do
	#Get basename for PBS
	fileID=$(echo $file | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	#set output file names
	MARKUP=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MARKUP.bam"}')
	METRICS=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_output_metrics.tab"}')
	FINAL=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_DEDUP.bam"}')

	sampleID=$fileID\_dedup.pbs
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output directory" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replciates" >> $sampleID
	echo "java -jar $PICARD MarkDuplicates  INPUT=$file OUTPUT=$MARKUP  METRICS_FILE=$METRICS  REMOVE_DUPLICATES='true' ASSUME_SORTED='true'  DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES'  READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100'  # Optional arguments  VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR" >> $sampleID
	echo "#remove duplicates, unmapped reads, etc." >> $sampleID
	echo "samtools view -bq 5 $MARKUP > $FINAL" >> $sampleID
	echo "#remove intermediate files" >> $sampleID
	echo "rm $MARKUP" >> $sampleID
	echo "#finish script" >> $sampleID

done
