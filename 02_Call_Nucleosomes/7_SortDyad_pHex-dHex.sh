#!/bin/bash

# Calculate pHex/dHex ratio across full-dyad +1 Nucleosomes and sort. Slice Top/Bottom 1K RefPTs.
# see 04_Figures/Fig_5c_d.sh

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/02_Call_Nucleosomes
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
KREBS=$WRK/../data/RefPT-Krebs/
GENOME=$WRK/../data/hg19_files/hg19.fa.fai

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

TEMP=Make_pHex-dHex
[ -d $TEMP ] || mkdir $TEMP

## =====Calculate pHex/dHex ratio=====

# Expand 2bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2 $KREBS/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad.bed -o $TEMP/PlusOne_2bp.bed

# Build pHex and dHex ranges
bedtools slop -i $TEMP/PlusOne_2bp.bed -g $GENOME -l 72 -r 18 -s > $TEMP/pHex.bed
bedtools slop -i $TEMP/PlusOne_2bp.bed -g $GENOME -l 18 -r 72 -s > $TEMP/dHex.bed

# Tag Pileup around Hex ranges for two insert size ranges
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -x 127 -n 92 --combined --cpu $THREADS $TEMP/pHex.bed $BAMFILE -M $TEMP/pHex
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -x 127 -n 92 --combined --cpu $THREADS $TEMP/dHex.bed $BAMFILE -M $TEMP/dHex

# Row-wise tally pileup signal
java -jar $SCRIPTMANAGER read-analysis aggregate-data -m -o $TEMP/pHex-dHex_combined.tab $TEMP/pHex_combined.cdt $TEMP/dHex_combined.cdt

# Calculate ratio:
#  - confirm IDs match up from paste
#  - exclude sites with 0 tags in either proximal or distal regions
#  - shuffle then resort by quotient (avoids subsort effects)
paste $KREBS/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad.bed <(sed '1d' $TEMP/pHex-dHex_combined.tab) \
	| awk '{OFS="\t"}{if ($4==$7 && $8!=0 && $9!=0) print $1,$2,$3,$4,($8/$9),$6}' \
	| shuf | sort -rnk5,5 \
	> $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad.bed

## =====Slice Info TSV into final BED files=====

# Slice out top and bottom 1k
head -n 1000 $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad.bed > $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-TOP-1K.bed
tail -n 1000 $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad.bed > $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-BOTTOM-1K.bed

[ -d $KREBS/200bp ] || mkdir $KREBS/200bp

# Expand 200bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad.bed -o $KREBS/200bp/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_200bp.bed

[ -d $KREBS/500bp ] || mkdir $KREBS/500bp

# Expand 500bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-TOP-1K.bed -o $KREBS/500bp/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-TOP-1K_500bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 500 $KREBS/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-BOTTOM-1K.bed -o $KREBS/500bp/PlusOneDyad_SORT-pHex-dHex_GROUP-Nuc-Dyad_GROUP-BOTTOM-1K_500bp.bed
