#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/Benzonase-ChIP_download-merge-dedup-merge.log.out
#SBATCH -e logs/Benzonase-ChIP_download-merge-dedup-merge.log.err

# Benzonase-ChIP and ChIP-exo data: script to hardcode the merging/renaming of PEGR BAM files into a standard file naming system

# data/BAM
#   |--ChIP-exo_Pol2_merge_hg38.bam
#   |--BNase-ChIP_H2A_merge_hg38.bam
#   |--BNase-ChIP_H2AZ_merge_hg38.bam
#   |--BNase-ChIP_H2B_merge_hg38.bam
#   |--BNase-ChIP_H3_merge_hg38.bam
#   |--BNase-ChIP_H3K27ac_merge_hg38.bam
#   |--BNase-ChIP_H3K27me3_merge_hg38.bam
#   |--BNase-ChIP_H3K36me3_merge_hg38.bam
#   |--BNase-ChIP_H3K4me1_merge_hg38.bam
#   |--BNase-ChIP_H3K4me3_merge_hg38.bam
#   |--BNase-ChIP_H3K9ac_merge_hg38.bam
#   |--BNase-ChIP_H4_merge_hg38.bam

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

[ -d ChIP ] || mkdir ChIP

cd ChIP

# Merge resequencing replicates
java -jar $PICARD MarkDuplicates -O H2A.Z_1.bam -M H2A.Z_1.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28453_H2A.Z_ab4174_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28794_H2A.Z_ab4174_K562_-_-_50Unuclease10min-0cycSonic_BX.bam
java -jar $PICARD MarkDuplicates -O H2A.Z_2.bam -M H2A.Z_2.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28461_H2A.Z_ab4174_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28799_H2A.Z_ab4174_K562_-_-_50Unuclease10min-0cycSonic_BX.bam
java -jar $PICARD MarkDuplicates -O H2B_1.bam -M H2B_1.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28454_H2B_ab52985_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28795_H2B_ab52985_K562_-_-_50Unuclease10min-0cycSonic_BX.bam
java -jar $PICARD MarkDuplicates -O H2B_2.bam -M H2B_2.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28462_H2B_ab52985_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28800_H2B_ab52985_K562_-_-_50Unuclease10min-0cycSonic_BX.bam
java -jar $PICARD MarkDuplicates -O H3K27me3_1.bam -M H3K27me3_1.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28456_H3K27me3_07-449_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28796_H3K27me3_07-449_K562_-_-_50Unuclease10min-0cycSonic_BX.bam
java -jar $PICARD MarkDuplicates -O H3K27me3_2.bam -M H3K27me3_2.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28464_H3K27me3_07-449_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28801_H3K27me3_07-449_K562_-_-_50Unuclease10min-0cycSonic_BX.bam
java -jar $PICARD MarkDuplicates -O H4_1.bam -M H4_1.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28457_H4_ab10158_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28797_H4_ab10158_K562_-_-_50Unuclease10min-0cycSonic_BX.bam
java -jar $PICARD MarkDuplicates -O H4_2.bam -M H4_2.metrics \
	--REMOVE_DUPLICATES 'true' --ASSUME_SORTED 'true' --DUPLICATE_SCORING_STRATEGY 'SUM_OF_BASE_QUALITIES' \
	--READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' --OPTICAL_DUPLICATE_PIXEL_DISTANCE '100' \
	-I 28465_H4_ab10158_K562_-_-_50Unuclease10min-0cycSonic_BX.bam \
	-I 28802_H4_ab10158_K562_-_-_50Unuclease10min-0cycSonic_BX.bam

# Rename non-resequenced replicates to standard naming
mv 19354_PolII_05-623B_K562_-_completeIMDM-37oC_-_XO.bam PolII_1.bam
mv 19355_PolII_ab76123_K562_-_completeIMDM-37oC_-_XO.bam PolII_2.bam

mv 28452_H2A_ab18255_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H2A_1.bam
mv 28460_H2A_ab18255_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H2A_2.bam

mv 25860_H3_ab1791_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3_1.bam
mv 25868_H3_ab1791_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3_2.bam
mv 25964_H3_ab1791_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3_3.bam
mv 25971_H3_ab1791_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3_4.bam
mv 28804_H3_ab1791_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3_5.bam
mv 28808_H3_ab1791_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3_6.bam

mv 25858_H3K27ac_ab4729_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K27ac_1.bam
mv 25866_H3K27ac_ab4729_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K27ac_2.bam
mv 25962_H3K27ac_ab4729_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K27ac_3.bam
mv 25969_H3K27ac_ab4729_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K27ac_4.bam
mv 28806_H3K27ac_ab4729_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K27ac_5.bam
mv 28810_H3K27ac_ab4729_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K27ac_6.bam

mv 25863_H3K36me3_ab9050_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K36me3_1.bam
mv 25871_H3K36me3_ab9050_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K36me3_2.bam
mv 25967_H3K36me3_ab9050_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K36me3_3.bam
mv 25974_H3K36me3_ab9050_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K36me3_4.bam

mv 28455_H3K4me1_ab8895_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me1_1.bam
mv 28463_H3K4me1_ab8895_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me1_2.bam

mv 25861_H3K4me3_ab8580_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me3_1.bam
mv 25869_H3K4me3_ab8580_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me3_2.bam
mv 25965_H3K4me3_ab8580_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me3_3.bam
mv 25972_H3K4me3_ab8580_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me3_4.bam
mv 28805_H3K4me3_ab8580_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me3_5.bam
mv 28809_H3K4me3_ab8580_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K4me3_6.bam

mv 25862_H3K9ac_07-352_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K9ac_1.bam
mv 25870_H3K9ac_07-352_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K9ac_2.bam
mv 25966_H3K9ac_07-352_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K9ac_3.bam
mv 25973_H3K9ac_07-352_K562_-_-_50Unuclease10min-0cycSonic_BX.bam H3K9ac_4.bam

#20min

# Pol2
java -jar $PICARD MergeSamFiles -O $BAMDIR/ChIP-exo_Pol2_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I PolII_1.bam -I PolII_2.bam
# H2A.Z
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H2AZ_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H2A.Z_1.bam -I H2A.Z_2.bam
# H2A
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H2A_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H2A_1.bam -I H2A_2.bam
# H2B
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H2B_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H2B_1.bam -I H2B_2.bam
# H3
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H3_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H3_1.bam -I H3_2.bam -I H3_3.bam \
	-I H3_4.bam -I H3_5.bam -I H3_6.bam
# H3K27ac
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H3K27ac_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H3K27ac_1.bam -I H3K27ac_2.bam -I H3K27ac_3.bam \
	-I H3K27ac_4.bam -I H3K27ac_5.bam -I H3K27ac_6.bam
# H3K27me3
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H3K27me3_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H3K27me3_1.bam -I H3K27me3_2.bam
# H3K36me3
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H3K36me3_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H3K36me3_1.bam -I H3K36me3_2.bam \
	-I H3K36me3_3.bam -I H3K36me3_4.bam
# H3K4me1
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H3K4me1_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H3K4me1_1.bam -I H3K4me1_2.bam
# H3K4me3
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H3K4me3_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H3K4me3_1.bam -I H3K4me3_2.bam -I H3K4me3_3.bam \
	-I H3K4me3_4.bam -I H3K4me3_5.bam -I H3K4me3_6.bam
# H3K9ac
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H3K9ac_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H3K9ac_1.bam -I H3K9ac_2.bam \
	-I H3K9ac_3.bam -I H3K9ac_4.bam
# H4
java -jar $PICARD MergeSamFiles -O $BAMDIR/BNase-ChIP_H4_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I H4_1.bam -I H4_2.bam

# Index
for FILE in $BAMDIR/BNase-ChIP*.bam;
do
	samtools index $FILE
done