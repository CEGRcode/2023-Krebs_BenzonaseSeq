# purpose - merge technical teplicates of benzonase-seq
# usage
# qq
#
# example
#
# 'qq'

#input directory
INPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/fastq_files

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files

#set hg19.fa
hg19=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/00_REFERENCE_DATA/GENOMES/hg19/hg19.fa

#set R1 and R2
R1=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/fastq_files/SRR16815400_1.fastq
R2=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/fastq_files/SRR16815400_2.fastq

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
#PBS -l walltime=10:00:00
#PBS -A open
cd $OUTPUT

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

	#Get basename for PBS
	R1a=$(echo $R1 | rev | cut -d"/" -f1 | rev | awk -F_ '{print $1}')
	fileID=$(echo $R1a | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
	echo $fileID
	#set output file names
	ALIGNED=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_ALIGNED.bam"}')
	MARKUP=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MARKUP.bam"}')
	METRICS=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_output_metrics.tab"}')
	FINAL=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_DEDUP.bam"}')
	MASTER=$(echo $fileID | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_master.bam"}')

	sampleID=$fileID\_align_dedup_filter.pbs
        rm -f $sampleID
	echo "$JOBSTATS" >> $sampleID
	echo "#set output directory" >> $sampleID
	echo "cd $OUTPUT" >> $sampleID
	echo "#align fastq files" >> $sampleID
	echo "bwa mem -t 4 -v 1 -T '30' -h '5' -M $hg19 $R1 $R2 | samtools sort -@4 -O bam -o $ALIGNED" >> $sampleID
	echo "#mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replciates" >> $sampleID
	echo "java -jar $PICARD MarkDuplicates  INPUT=$ALIGNED OUTPUT=$MARKUP  METRICS_FILE=$METRICS  REMOVE_DUPLICATES='true' ASSUME_SORTED='true'  DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES'  READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100'  # Optional arguments  VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR" >> $sampleID
	echo "#remove duplicates, unmapped reads, etc." >> $sampleID
	echo "samtools view -bq 5 $MARKUP > $FINAL" >> $sampleID
	echo "#filter resulting bam to match PEGR" >> $sampleID
	echo "samtools view -o $MASTER -h -b -f 0x1 -F 0x404 $FINAL" >> $sampleID
	echo "#remove intermediate file" >> $sampleID
	echo "rm $ALIGNED" >> $sampleID
	echo "rm $MARKUP" >> $sampleID
	echo "rm $FINAL" >> $sampleID
	echo "#finish script" >> $sampleID
