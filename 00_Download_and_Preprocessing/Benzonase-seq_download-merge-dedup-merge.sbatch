#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 6:00:00
#SBATCH -A open
#SBATCH -o logs/Benzonase-seq_download-merge-dedup-merge.log.out
#SBATCH -e logs/Benzonase-seq_download-merge-dedup-merge.log.err

# Script to hardcode the merging/renaming of PEGR BAM files into a standard file naming system

# data/BAM
#   |--BNase-seq_50U-3min_merge_hg38.bam
#   |--BNase-seq_50U-30min_merge_hg38.bam
#   |--BNase-seq_50U-10min_merge_hg38.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Script shortcuts
PICARD=$WRK/../bin/picard.jar
BAMDIR=$WRK/../data/BAM

[ -d Benzonase ] || mkdir Benzonase

cd Benzonase

# Merge and MarkDuplicates on resequencing replicates (biorep 12-17)
# Mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replicates
# Optional arguments: VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR
java -jar $PICARD MarkDuplicates -O 50U-10min_12.bam -M 50U-10min_12.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33023_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam \
	-I 33159_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam
java -jar $PICARD MarkDuplicates -O 50U-10min_13.bam -M 50U-10min_13.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33024_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam \
	-I 33160_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-total-input_BI.bam
java -jar $PICARD MarkDuplicates -O 50U-10min_14.bam -M 50U-10min_14.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33025_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam \
	-I 33161_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam
java -jar $PICARD MarkDuplicates -O 50U-10min_15.bam -M 50U-10min_15.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33026_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam \
	-I 33162_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-input_BI.bam
java -jar $PICARD MarkDuplicates -O 50U-10min_16.bam -M 50U-10min_16.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33027_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam \
	-I 33163_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam
java -jar $PICARD MarkDuplicates -O 50U-10min_17.bam -M 50U-10min_17.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 33028_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam \
	-I 33164_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-pellet-input_BI.bam

# Rename non-resequenced replicates to standard naming
mv 28385_Input_-_K562_-_-_50Unuclease3min-0cycSonic_BI.bam 50U-3min_1.bam
mv 28386_Input_-_K562_-_-_50Unuclease3min-0cycSonic_BI.bam 50U-3min_2.bam
mv 28393_Input_-_K562_-_-_50Unuclease30min-0cycSonic_BI.bam 50U-30min_1.bam
mv 28394_Input_-_K562_-_-_50Unuclease30min-0cycSonic_BI.bam 50U-30min_2.bam
mv 25769_Input_-_K562_-_-_50Unuclease-0cycSonic_BI.bam 50U-10min_1.bam
mv 25976_Input_-_K562_-_-_50Unuclease10min-0cycSonic_BI.bam 50U-10min_2.bam
mv 28326_Input_-_K562_-_-_50Unuclease10min-0cycSonic-newQuenchingbuffer_BI.bam 50U-10min_3.bam
mv 28389_Input_-_K562_-_-_50Unuclease10min-0cycSonic_BI.bam 50U-10min_4.bam
mv 28390_Input_-_K562_-_-_50Unuclease10min-0cycSonic_BI.bam 50U-10min_5.bam
mv 28468_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Total_BI.bam 50U-10min_6.bam
mv 28469_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Pellet_BI.bam 50U-10min_7.bam
mv 28471_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam 50U-10min_8.bam
mv 28472_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam 50U-10min_9.bam
mv 28473_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam 50U-10min_10.bam
mv 28820_Input_-_K562_-_-_50Unuclease10min-0cycSonic-Super_BI.bam 50U-10min_11.bam
mv 33646_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-EndRepair-Atail-ATligation_BI.bam 50U-10min_18.bam
mv 33647_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-noEndRepair-Atail-ATligation_BI.bam 50U-10min_19.bam
mv 33648_Input_-_K562_-_-_XL-50Unuclease-0cycSonic-super-EndRepair-Atail-ATligation_BI.bam 50U-10min_20.bam

# 50U - 3min
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-seq_50U-3min_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I 50U-3min_1.bam -I 50U-3min_2.bam
# 50U - 30min
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-seq_50U-30min_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I 50U-30min_1.bam -I 50U-30min_2.bam
# 50U - 10min
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I 50U-10min_1.bam -I 50U-10min_2.bam -I 50U-10min_3.bam \
	-I 50U-10min_4.bam -I 50U-10min_5.bam -I 50U-10min_6.bam \
	-I 50U-10min_7.bam -I 50U-10min_8.bam -I 50U-10min_9.bam \
	-I 50U-10min_10.bam -I 50U-10min_11.bam -I 50U-10min_12.bam \
	-I 50U-10min_13.bam -I 50U-10min_14.bam -I 50U-10min_15.bam \
	-I 50U-10min_16.bam -I 50U-10min_17.bam -I 50U-10min_18.bam \
	-I 50U-10min_19.bam -I 50U-10min_20.bam

samtools index $BAMDIR/BNase-seq_50U-3min_merge_hg38.bam
samtools index $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam
samtools index $BAMDIR/BNase-seq_50U-30min_merge_hg38.bam
