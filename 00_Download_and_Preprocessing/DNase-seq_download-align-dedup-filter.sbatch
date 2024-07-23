#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 3:00:00
#SBATCH -A open
#SBATCH -o logs/1_download_align_dedup_DNase.log.out
#SBATCH -e logs/1_download_align_dedup_DNase.log.err

# Download, align, and dedup the larger replicate of DNase-seq.
# see /storage/group/bfp2/default/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/jobs/b_align_dedup_filter_230815.sh

# data/BAM
#   |--DNase-seq_-_rep1_hg19.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
###

# Dependencies
# - bwa
# - java
# - samtools
# - SRA toolkit (fasterq-dump)
# - wget

set -exo
module load samtools
cd $WRK

# Inputs and outputs
BAMDIR=$WRK/../data/BAM
GENOME=$WRK/../data/hg19_files/hg19.fa

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d FQ ] || mkdir FQ
[ -d BAM ] || mkdir BAM

SRR=SRR16815400

# Download raw FASTQ file
fasterq-dump -O FQ --split-files $SRR

# Align each sample
bwa mem -t 4 -v 1 -T '30' -h '5' -M $GENOME FQ/$SRR\_1.fastq FQ/$SRR\_2.fastq | samtools sort - -@4 -o BAM/$SRR\_sorted.bam

# Mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replciates
# Optional arguments: VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR
java -jar $PICARD MarkDuplicates INPUT=BAM/$SRR\_sorted.bam OUTPUT=BAM/$SRR\_markdup.bam METRICS_FILE=BAM/$SRR\_metrics.txt \
	REMOVE_DUPLICATES='true' ASSUME_SORTED='true' DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' \
	READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100'

# Remove duplicates, unmapped reads, with mapping quality threshold
samtools view -bq 5 BAM/$SRR\_markdup.bam -o BAM/$SRR\_dedup.bam

# Filter resulting BAM
samtools view -h -b -f 0x1 -F 0x404 BAM/$SRR\_dedup.bam -o $BAMDIR/DNase-seq_-_rep1_hg19.bam
