#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/MNase-seq-Titrations_download-align-dedup-filter.log.out
#SBATCH -e logs/MNase-seq-Titrations_download-align-dedup-filter.log.err

# Download, align, and dedup MNase sequencing data from an MNase titration.
# see /storage/group/bfp2/default/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/jobs/b_align_dedup_filter_230815.sh

# data/BAM
#   |--MNase-seq_5U_rep1_hg19.bam
#   |--MNase-seq_21U_rep1_hg19.bam
#   |--MNase-seq_79U_rep1_hg19.bam
#   |--MNase-seq_304U_rep1_hg19.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
THREADS=16
###

# Dependencies
# - bwa
# - java
# - samtools
# - SRA toolkit (fasterq-dump)
# - wget

set -exo
module load bwa
module load samtools
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/align

# Inputs and outputs
BAMDIR=$WRK/../data/BAM
GENOME=$WRK/../data/hg19_files/hg19.fa

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d FQ ] || mkdir FQ
[ -d BAM ] || mkdir BAM

for SRR in "SRR3211679" "SRR3211680" "SRR3211681" "SRR3211682";
do
	# Download raw FASTQ file
	fasterq-dump -O FQ --split-files $SRR

	# Compress FASTQ files
	pigz FQ/${SRR}_1.fastq
	pigz FQ/${SRR}_2.fastq

	# Align each sample
	bwa mem -t $THREADS -v 1 $GENOME FQ/${SRR}_1.fastq FQ/${SRR}_2.fastq | samtools sort -@ $THREADS > BAM/${SRR}_sorted.bam

	# Mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replicates
	# Optional arguments: VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR
	java -jar $PICARD MarkDuplicates -I BAM/${SRR}_sorted.bam -O BAM/${SRR}_dedup.bam -M BAM/${SRR}_dedup.metrics \
		--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
		--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100'
		--VALIDATION_STRINGENCY 'LENIENT' --TAGGING_POLICY All --QUIET true --VERBOSITY ERROR
done

# Filter resulting BAM to match PEGR & format final BAM filename
samtools view -b -f 0x1 -F 0x404 BAM/SRR3211679_dedup.bam -o $BAMDIR/MNase-seq_5U_rep1_hg19.bam
samtools view -b -f 0x1 -F 0x404 BAM/SRR3211680_dedup.bam -o $BAMDIR/MNase-seq_21U_rep1_hg19.bam
samtools view -b -f 0x1 -F 0x404 BAM/SRR3211681_dedup.bam -o $BAMDIR/MNase-seq_79U_rep1_hg19.bam
samtools view -b -f 0x1 -F 0x404 BAM/SRR3211682_dedup.bam -o $BAMDIR/MNase-seq_304U_rep1_hg19.bam
