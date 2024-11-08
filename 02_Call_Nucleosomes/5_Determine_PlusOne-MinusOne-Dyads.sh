#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/5_Determine_PlusOne-MinusOne-Dyads.log.out
#SBATCH -e logs/5_Determine_PlusOne-MinusOne-Dyads.log.err

# Determine closest inferred nucleosome to CoPRO-determined TSS
# see 03_NucCalls/01a_NucCall_expressed_TSS.sh
# see 03_NucCalls/01b_NucCall_unexpressed_TSS.sh
# see 03_NucCalls/02_NFR.sh

# data/RefPT-Krebs
#   |--MinusOneDyad_SORT-DistToExpressedTSS.bed
#   |--PlusOneDyad_SORT-DistToExpressedTSS.bed
#   |--MinusOneDyad_SORT-DistToUnexpressedTSS.bed
#   |--PlusOneDyad_SORT-DistToUnexpressedTSS.bed
#   |--MinusOneDyad_SORT-Expression.bed
#   |--PlusOneDyad_SORT-Expression.bed
#   |--PlusOneDyad_SORT-Expression_WithUnexpressed.bed
#   |--PlusOneDyad_SORT-DistToExpressedTSS_GROUP-CpGIsland-NoOverlap.bed
#   |--PlusOneDyad_SORT-DistToExpressedTSS_GROUP-CpGIsland-Overlap.bed
#   |--2000bp
#     |--PlusOneDyad_SORT-Expression_2000bp.bed
#     |--PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp.bed


### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/02_Call_Nucleosomes
###

# Dependencies
# - java
# - perl
# - python

set -exo
module load anaconda3
module load bedtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
OTHER=$WRK/../data/RefPT-Other/
KREBS=$WRK/../data/RefPT-Krebs/
NUCLEOSOME=$KREBS/BNase-Nucleosomes.bed
GENOME=$WRK/../data/hg38_files/hg38.fa.fai
CPG=../data/RefPT-Other/CpGIslands.bed
EXPRESSED=$KREBS/TSS_GROUP-Expressed_SORT-Expression.bed
UNEXPRESSED=$KREBS/TSS_GROUP-Unexpressed.bed

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
FILTERL=../bin/filter_BED_by_list_ColumnSelect.pl

TEMP=MakePlusMinus
[ -d $TEMP ] || mkdir $TEMP

# Update score column with rank RNA expression/"unexpression" (for resorting later). Assumes starting sort is RNA expresssion.
awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,NR,$6}' $EXPRESSED > $TEMP/TSS_SORT-RankExpression.bed
awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,$4,NR,$6}' $UNEXPRESSED > $TEMP/uTSS_SORT-RankSort.bed

# Expand 2bp (expressed/unexpressed)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $TEMP/TSS_SORT-RankExpression.bed -o $TEMP/TSS_SORT-RankExpression_1bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1 $TEMP/uTSS_SORT-RankSort.bed -o $TEMP/uTSS_SORT-RankSort_1bp.bed

# Genomic sort
bedtools sort -i $NUCLEOSOME > $TEMP/Nucleosomes.bed
bedtools sort -i $TEMP/TSS_SORT-RankExpression_1bp.bed > $TEMP/TSS_SORT-Genomic.bed
bedtools sort -i $TEMP/uTSS_SORT-RankSort_1bp.bed > $TEMP/uTSS_SORT-Genomic.bed

# Get closest up/downstream inferred nucleosomes (get multiple closest hits for minus one)
bedtools closest -k 4 -t first -D a -d -id -a $TEMP/TSS_SORT-Genomic.bed  -b $TEMP/Nucleosomes.bed > $TEMP/TSS_upstream_Nucleosomes.bed
bedtools closest -k 1 -t first -D a -d -iu -a $TEMP/TSS_SORT-Genomic.bed  -b $TEMP/Nucleosomes.bed > $TEMP/TSS_downstream_Nucleosomes.bed
bedtools closest -k 4 -t first -D a -d -id -a $TEMP/uTSS_SORT-Genomic.bed -b $TEMP/Nucleosomes.bed > $TEMP/uTSS_upstream_Nucleosomes.bed
bedtools closest -k 1 -t first -D a -d -iu -a $TEMP/uTSS_SORT-Genomic.bed -b $TEMP/Nucleosomes.bed > $TEMP/uTSS_downstream_Nucleosomes.bed

# (Upstream) Sort by distance (TSS to Nuc) - Filter to exclude distances farther than 5kb - Filter to keep closest non-overlapping Dyad for MinusOne
sort -rnk13,13 $TEMP/TSS_upstream_Nucleosomes.bed \
	| awk '{FS="\t"}{OFS="\t"}{if ($13!=0 && $13>-5000 && $7!=".") print;}' \
	| sort -uk5,5 | sort -rnk13,13 \
	> $TEMP/MinusOneDyad_SORT-DistToExpressedTSS.tsv
sort -rnk13,13 $TEMP/uTSS_upstream_Nucleosomes.bed \
	| awk '{OFS="\t"}{FS="\t"}{if ($13!=0 && $13>-5000 && $7!=".") print;}' \
	| sort -uk5,5 | sort -rnk13,13 \
	> $TEMP/MinusOneDyad_SORT-DistToUnexpressedTSS.tsv

# (Downstream) Sort by distance (TSS to Nuc) - Filter to exclude distances farther than 5kb - Can be overlapping Dyad for PlusOne
sort -nk13,13 $TEMP/TSS_downstream_Nucleosomes.bed \
	| awk '{OFS="\t"}{FS="\t"}{if ($13<5000 && $7!=".") print;}' \
	> $TEMP/PlusOneDyad_SORT-DistToExpressedTSS.tsv
sort -nk13,13 $TEMP/uTSS_downstream_Nucleosomes.bed \
	| awk '{OFS="\t"}{FS="\t"}{if ($13<5000 && $7!=".") print;}' \
	> $TEMP/PlusOneDyad_SORT-DistToUnexpressedTSS.tsv


## =====Slice Info TSV into BED=====

# Slice columns for proper Nuc BED (id from Nuc, strand from TSS, score from dist)
awk '{OFS="\t"}{FS="\t"}{print $7,$8,$9,$10,$13,$6}' $TEMP/MinusOneDyad_SORT-DistToExpressedTSS.tsv   > $KREBS/MinusOneDyad_SORT-DistToExpressedTSS.bed
awk '{OFS="\t"}{FS="\t"}{print $7,$8,$9,$10,$13,$6}' $TEMP/PlusOneDyad_SORT-DistToExpressedTSS.tsv    > $KREBS/PlusOneDyad_SORT-DistToExpressedTSS.bed
awk '{OFS="\t"}{FS="\t"}{print $7,$8,$9,$10,$13,$6}' $TEMP/MinusOneDyad_SORT-DistToUnexpressedTSS.tsv > $KREBS/MinusOneDyad_SORT-DistToUnexpressedTSS.bed
awk '{OFS="\t"}{FS="\t"}{print $7,$8,$9,$10,$13,$6}' $TEMP/PlusOneDyad_SORT-DistToUnexpressedTSS.tsv  > $KREBS/PlusOneDyad_SORT-DistToUnexpressedTSS.bed

# Sort by TSS expresssion for proper Nuc BED (id from Nuc, strand from TSS, score from TSS)
sort -nk5,5 $TEMP/MinusOneDyad_SORT-DistToExpressedTSS.tsv \
	| awk '{OFS="\t"}{FS="\t"}{print $7,$8,$9,$10,$5,$6}' \
	> $KREBS/MinusOneDyad_SORT-Expression.bed
sort -nk5,5 $TEMP/PlusOneDyad_SORT-DistToExpressedTSS.tsv \
	| awk '{OFS="\t"}{FS="\t"}{print $7,$8,$9,$10,$5,$6}' \
	> $KREBS/PlusOneDyad_SORT-Expression.bed

# Subset by Dyad calls derived from full-length octomer fragments
cut -f4 $KREBS/PlusOneDyad_SORT-Expression.bed | grep '^Nuc-' > $TEMP/Nuc-Dyad.ids
perl $FILTERL $KREBS/PlusOneDyad_SORT-Expression.bed $TEMP/Nuc-Dyad.ids 3 keep $KREBS/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad.bed

# PlusOneDyad sorted by RankExpression with unexpressed PlusOneDyad
cat $KREBS/PlusOneDyad_SORT-Expression.bed $KREBS/PlusOneDyad_SORT-DistToUnexpressedTSS.bed > $KREBS/PlusOneDyad_SORT-Expression_WithUnexpressed.bed

# Expand 2000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/PlusOneDyad_SORT-DistToExpressedTSS.bed -o $KREBS/2000bp/PlusOneDyad_SORT-DistToExpressedTSS_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/MinusOneDyad_SORT-DistToExpressedTSS.bed -o $KREBS/2000bp/MinusOneDyad_SORT-DistToExpressedTSS_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/PlusOneDyad_SORT-DistToUnexpressedTSS.bed -o $KREBS/2000bp/PlusOneDyad_SORT-DistToUnexpressedTSS_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/MinusOneDyad_SORT-DistToUnexpressedTSS.bed -o $KREBS/2000bp/MinusOneDyad_SORT-DistToUnexpressedTSS_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/PlusOneDyad_SORT-Expression.bed -o $KREBS/2000bp/PlusOneDyad_SORT-Expression_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/MinusOneDyad_SORT-Expression.bed -o $KREBS/2000bp/MinusOneDyad_SORT-Expression_2000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/PlusOneDyad_SORT-Expression_WithUnexpressed.bed -o $KREBS/2000bp/PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp.bed


## =====Group PlusOne Dyad by CpG island overlap=====

# Intersect CpG Islands with PlusOne
bedtools intersect    -wa -a $KREBS/PlusOneDyad_SORT-DistToExpressedTSS.bed -b $CPG > $KREBS/PlusOneDyad_SORT-DistToExpressedTSS_GROUP-CpGIsland-Overlap.bed
bedtools intersect -v -wa -a $KREBS/PlusOneDyad_SORT-DistToExpressedTSS.bed -b $CPG > $KREBS/PlusOneDyad_SORT-DistToExpressedTSS_GROUP-CpGIsland-NoOverlap.bed

# Sort by CpGIsland length and dedup PlusOne
sort -rnk11,11 $TEMP/PlusOne-CpGIsland_Intersect.tsv | sort -uk4,4 | sort -rnk11,11 > $KREBS/TSS_GROUP-Expressed_SORT-CpG.bed

# Expand 2000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 2000 $KREBS/TSS_GROUP-Expressed_SORT-CpG.bed -o $KREBS/2000bp/TSS_GROUP-Expressed_SORT-CpG_2000bp.bed