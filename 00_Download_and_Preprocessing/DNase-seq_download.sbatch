#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=4gb
#SBATCH -t 0:10:00
#SBATCH -A open
#SBATCH -o logs/DNase-seq_download.log.out
#SBATCH -e logs/DNase-seq_download.log.err

# Download and rename the ENCODE DNase-seq BAM file
# see /storage/group/bfp2/default/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/jobs/b_align_dedup_filter_230815.sh

# data/BAM
#   |--DNase-seq_ENCFF425WDA_rep1_hg38.bam
#   |--DNase-seq_ENCFF518XTC_rep1_hg38.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
###

# Dependencies
# - samtools
# - wget

set -exo
module load samtools

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Low depth (<20 mil) from ENCSR000EOT
wget -c -O $BAMDIR/DNase-seq_ENCFF425WDA_rep1_hg38.bam https://www.encodeproject.org/files/ENCFF425WDA/@@download/ENCFF425WDA.bam
# GM with TALEN (DHS sites)
wget -c -O $BAMDIR/DNase-seq_ENCFF518XTC_rep1_hg38.bam https://www.encodeproject.org/files/ENCFF518XTC/@@download/ENCFF518XTC.bam

samtools index $BAMDIR/DNase-seq_ENCFF425WDA_rep1_hg38.bam
samtools index $BAMDIR/DNase-seq_ENCFF518XTC_rep1_hg38.bam