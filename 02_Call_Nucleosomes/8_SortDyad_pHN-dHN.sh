#!/bin/bash

# Calculate pHN/dHN ratio across full-dyad +1 Nucleosomes and sort. Slice Top/Bottom 2,500 RefPTs (skip over very bottom ~500 sites).
# see 04_Figures/Fig_5c_d.sh

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
THREADS=4
###

# Dependencies
# - java
# - python

set -exo
module load anaconda3
module load bedtools
source activate bx

# Inputs and outputs
BAMFILE=$WRK/../data/BAM/BNase-ChIP_H3K4me3_merge_hg19.bam
KREBS=$WRK/../data/RefPT-Krebs
GENOME=$WRK/../data/hg19_files/hg19.fa.fai

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

TEMP=Make_pHN-dHN
[ -d $TEMP ] || mkdir $TEMP

## =====Calculate pHN/dHN ratio=====

# Expand 2bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2 $KREBS/PlusOneDyad_SORT-Expression.bed -o $TEMP/PlusOne_2bp.bed

# Build pHN and dHN ranges
bedtools slop -i $TEMP/PlusOne_2bp.bed -g $GENOME -l 82 -r -4 -s > $TEMP/pHN.bed
bedtools slop -i $TEMP/PlusOne_2bp.bed -g $GENOME -l -6 -r 84 -s > $TEMP/dHN.bed

# Tag Pileup around HN ranges for two insert size ranges
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -x 80 --combined --cpu $THREADS $TEMP/pHN.bed $BAMFILE -M $TEMP/pHN
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -x 80 --combined --cpu $THREADS $TEMP/dHN.bed $BAMFILE -M $TEMP/dHN

# Row-wise tally pileup signal
java -jar $SCRIPTMANAGER read-analysis aggregate-data -m -o $TEMP/pHN-dHN_combined.tab $TEMP/pHN_combined.cdt $TEMP/dHN_combined.cdt

# Calculate ratio:
#  - confirm IDs match up from paste
#  - exclude sites with 0 tags in either proximal or distal regions
#  - shuffle then resort by quotient (avoids subsort effects)
paste $KREBS/PlusOneDyad_SORT-Expression.bed <(sed '1d' $TEMP/pHN-dHN_combined.tab) \
	| awk '{OFS="\t"}{if ($4==$7 && $8!=0 && $9!=0) print $1,$2,$3,$4,($8/$9),$6}' \
	| shuf | sort -rnk5,5 \
	> $KREBS/PlusOneDyad_SORT-pHN-dHN.bed

## =====Slice Info TSV into final BED files=====

# Slice out top and bottom 2.5k (bottom doesn't include very last 573 sites)
head -n 2500 $KREBS/PlusOneDyad_SORT-pHN-dHN.bed > $KREBS/PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500.bed
head -n -573 $KREBS/PlusOneDyad_SORT-pHN-dHN.bed | tail -n 2500 > $KREBS/PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500.bed

[ -d $KREBS/2000bp ] || mkdir $KREBS/2000bp

# Expand 2000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/PlusOneDyad_SORT-pHN-dHN.bed -o $KREBS/2000bp/PlusOneDyad_SORT-pHN-dHN_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500.bed -o $KREBS/2000bp/PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500.bed -o $KREBS/2000bp/PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500_2000bp.bed
