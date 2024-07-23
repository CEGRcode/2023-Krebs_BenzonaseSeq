#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/MNase-seq-ENCODE_download_merge.log.out
#SBATCH -e logs/MNase-seq-ENCODE_download_merge.log.err

# Download deep MNase sequencing (single-end, ABI SOLiD) and merge from ENCODE
# No need to MarkDuplicates b/c ENCODE says these are bio replicates (and
# pipeline already filters)

# data/BAM
#   |--MNase-seq_ENCODE_merge_hg19.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
###

# Dependencies
# - java
# - samtools
# - wget

set -exo
module load bwa
module load samtools
cd $WRK

# Inputs and outputs
BAMDIR=$WRK/../data/BAM
GENOME=$WRK/../data/hg19_files/hg19.fa.colorspace

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d ENCODE ] || mkdir ENCODE

# Download filtered BAM files
#wget -c -O ENCODE/ENCFF000VMJ.bam "https://www.encodeproject.org/files/ENCFF000VMJ/@@download/ENCFF000VMJ.bam"
wget -c -O ENCODE/ENCFF000VNR.bam "https://www.encodeproject.org/files/ENCFF000VNR/@@download/ENCFF000VNR.bam"       #(?)
#wget -c -O ENCODE/ENCFF000VML.bam "https://www.encodeproject.org/files/ENCFF000VML/@@download/ENCFF000VML.bam"
wget -c -O ENCODE/ENCFF000VMM.bam "https://www.encodeproject.org/files/ENCFF000VMM/@@download/ENCFF000VMM.bam"
wget -c -O ENCODE/ENCFF000VMN.bam "https://www.encodeproject.org/files/ENCFF000VMN/@@download/ENCFF000VMN.bam"
wget -c -O ENCODE/ENCFF000VMO.bam "https://www.encodeproject.org/files/ENCFF000VMO/@@download/ENCFF000VMO.bam"
wget -c -O ENCODE/ENCFF000VMP.bam "https://www.encodeproject.org/files/ENCFF000VMP/@@download/ENCFF000VMP.bam"
wget -c -O ENCODE/ENCFF000VNQ.bam "https://www.encodeproject.org/files/ENCFF000VNQ/@@download/ENCFF000VNQ.bam"
wget -c -O ENCODE/ENCFF000VNP.bam "https://www.encodeproject.org/files/ENCFF000VNP/@@download/ENCFF000VNP.bam"
wget -c -O ENCODE/ENCFF000VNO.bam "https://www.encodeproject.org/files/ENCFF000VNO/@@download/ENCFF000VNO.bam"
wget -c -O ENCODE/ENCFF000VMK.bam "https://www.encodeproject.org/files/ENCFF000VMK/@@download/ENCFF000VMK.bam"

# Merge BAM files
samtools merge -f $BAMDIR/MNase-seq_ENCODE_merge_hg19.bam \
	ENCODE/ENCFF000VMM.bam ENCODE/ENCFF000VMN.bam \
	ENCODE/ENCFF000VMO.bam ENCODE/ENCFF000VMP.bam \
	ENCODE/ENCFF000VNQ.bam ENCODE/ENCFF000VNP.bam \
	ENCODE/ENCFF000VNO.bam ENCODE/ENCFF000VMK.bam

# Index
samtools index $BAMDIR/MNase-seq_ENCODE_merge_hg19.bam

# Clean-up
rm -r ENCODE
