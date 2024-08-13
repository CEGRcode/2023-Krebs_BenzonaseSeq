#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 10:00:00
#SBATCH -A open
#SBATCH -o logs/CUTRUN_download-align-dedup-merge.log.out-%a
#SBATCH -e logs/CUTRUN_download-align-dedup-merge.log.err-%a
#SBATCH --array 1-8

# SLURM_ARRAY_TASK_ID=2

# Download CUT&RUN raw data (paired-end), merge, trim, align, dedup, filter, merge

# data/BAM
#   |--CUTRUN_CTCF_merge_hg38.bam
#   |--CUTRUN_H2AZ_merge_hg38.bam
#   |--CUTRUN_H3K4me1_merge_hg38.bam
#   |--CUTRUN_H3K4me3_merge_hg38.bam
#   |--CUTRUN_H3K27ac_merge_hg38.bam
#   |--CUTRUN_H3K27ac-Millipore_merge_hg38.bam
#   |--CUTRUN_H3K27me3_merge_hg38.bam
#   |--CUTRUN_IgG_merge_hg38.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
THREADS=8
###

# Dependencies
# - SRA toolkit (fasterq-dump)
# - java
# - samtools

set -exo
module load anaconda3
module load samtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/align

# Downloaded from https://raw.githubusercontent.com/timflutre/trimmomatic
PEADAPT=Trimmomatic-0.36/adapters/TruSeq3-PE.fa
GENOME=../data/hg38_files/hg38.fa
BAMDIR=../data/BAM

PICARD=../bin/picard.jar
TRIM=Trimmomatic-0.36/trimmomatic-0.36.jar

# CUT&RUN from 4DNucleome
[ -d CUTRUN ] || mkdir CUTRUN

# Identify metadata file and create its own directory
METADATA=`ls CUTRUN/*.tsv | head -n ${SLURM_ARRAY_TASK_ID} | tail -1`
EXPERIMENT=`basename $METADATA ".tsv"`
ODIR=CUTRUN/$EXPERIMENT

# Parse target & (remove "." for H2A.Z)
TARGET=`echo $EXPERIMENT | awk -F"_" '{print $1}'`
TARGET=${TARGET/./}

[ -d $ODIR ] || mkdir $ODIR

# Parse out header and comments from metadata file
TEMP=$ODIR/temp.txt
grep -v '^[#(File)]' $METADATA > $TEMP

# Count N elements in metadata
NCOUNT=`wc -l $TEMP | awk '{print $1}'`

# Download  each FASTQ file (format filename using metadata)
for N in $(seq 1 $NCOUNT);
do
	URL=`head -n $N $TEMP | tail -1 | awk 'BEGIN{FS="\t"}{print $31}'`
	FDNFF=`basename $URL ".fastq.gz"`
	EXP=`head -n $N $TEMP | tail -1 | awk '{print $3}'`
	# TARGET=`head -n $N $TEMP | tail -1 | awk 'BEGIN{FS="\t"}{gsub(" protein", "", $21); print $21 }'`
	READ=`head -n $N $TEMP | tail -1 | awk 'BEGIN{FS="\t"}{print $15}'`

	# Set Read2 files to use Read1 file id
	[ $READ -eq 2 ] && FDNFF=`head -n $N $TEMP | tail -1 | awk 'BEGIN{FS="\t"}{print $14}'`

	# Download FASTQ file
	wget -O $ODIR/$FDNFF.$EXP\_R$READ.fastq.gz "${URL}"
done

# List out experiments
EXPLIST=`ls $ODIR/4DN*.fastq.gz |cut -d"." -f2 |cut -d"_" -f1 |sort |uniq`

# Merge, trim, align, dedup, filter
for EXP in ${EXPLIST[@]};
do
	BASE=`basename $EXP "_R1.fastq.gz"`
	BASE=$ODIR/Merged_$BASE

	# Merge sequencing replicates (FASTQ)
	cat $ODIR/*.$EXP\_R1.fastq.gz > $BASE\_R1.fastq.gz
	cat $ODIR/*.$EXP\_R2.fastq.gz > $BASE\_R2.fastq.gz

	# Trim each set of merged reads (FASTQ)
	java -jar $TRIM PE -threads $THREADS -phred33 $BASE\_R1.fastq.gz $BASE\_R2.fastq.gz $BASE\_trim_R1.fastq.gz $BASE\_fail_R1.fastq.gz $BASE\_trim_R2.fastq.gz $BASE\_fail_R2.fastq.gz ILLUMINACLIP:$PEADAPT:2:15:4:1:true MINLEN:20

	# Align and sort
	bowtie2 --dovetail --threads $THREADS -x $GENOME -1 $BASE\_trim_R1.fastq.gz -2 $BASE\_trim_R2.fastq.gz | samtools sort > $BASE.bam

	# Mark duplicates and filter
	java -jar $PICARD MarkDuplicates -I $BASE.bam -M $BASE\_dedup.metrics -O $BASE\_dedup.bam --VALIDATION_STRINGENCY LENIENT
	samtools view -F 1024 -f 2 -b $BASE\_dedup.bam > $BASE\_final.bam

	# Build picard merge args string
	MERGE_LIST="-I ${BASE}_final.bam $MERGE_LIST"
done

# Merge target reps
java -jar $PICARD MergeSamFiles -O $BAMDIR/CUTRUN_$TARGET\_merge_hg38.bam $MERGE_LIST
samtools index $BAMDIR/CUTRUN_$TARGET\_merge_hg38.bam
