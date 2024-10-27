#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/CoPRO_download-ENCODE.log.out
#SBATCH -e logs/CoPRO_download-ENCODE.log.err

# Download capped CoPRO data from ENCODE and merge

# data/BAM
#   |--MNase-seq_ENCODE_merge_hg38.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
THREADS=8
###

# Dependencies
# - java
# - samtools
# - wget

set -exo
module load samtools
cd $WRK

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d ENCODE ] || mkdir ENCODE

# Download BAM files
wget -c -O ENCODE/ENCFF663UAN.bam https://www.encodeproject.org/files/ENCFF663UAN/@@download/ENCFF663UAN.bam
wget -c -O ENCODE/ENCFF877SQU.bam https://www.encodeproject.org/files/ENCFF877SQU/@@download/ENCFF877SQU.bam

# Merge BAM files
java -jar $PICARD MergeSamFiles -O $BAMDIR/CoPRO_Capped_merge_hg38.bam \
	-I ENCODE/ENCFF663UAN.bam \
	-I ENCODE/ENCFF877SQU.bam

# Index
samtools index $BAMDIR/CoPRO_Capped_merge_hg38.bam

# Clean-up
rm -r ENCODE