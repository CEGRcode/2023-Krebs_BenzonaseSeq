#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/GROcap_download-ENCODE.log.out
#SBATCH -e logs/GROcap_download-ENCODE.log.err

# Download capped GRO-cap data from ENCODE

# data/BAM
#   |--GRO-cap_ENCODE_merge_hg38.bam

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

# Download BAM files
wget -c -O $BAMDIR/GRO-cap_ENCODE_merge_hg38.bam https://www.encodeproject.org/files/ENCFF112MXH/@@download/ENCFF112MXH.bam

# Index
samtools index $BAMDIR/GRO-cap_ENCODE_merge_hg38.bam