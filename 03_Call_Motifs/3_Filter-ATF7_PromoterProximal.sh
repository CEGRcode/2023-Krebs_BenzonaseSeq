#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:05:00
#SBATCH -A open
#SBATCH -o logs/3_Filter-ATF7_PromoterProximal.log.out
#SBATCH -e logs/3_Filter-ATF7_PromoterProximal.log.err

# Re-intersect all ATF7 motif with ChIP-seq peaks and sort into
# promoter-proximal and not promoter-proximal (NFR overlap or not)

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/03_Call_Motifs
###

# Dependencies
# - bedtools
# - java

set -exo
module load anaconda3
module load bedtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif
NFR_BEDFILE=$WRK/../data/RefPT-Krebs/NFR_SORT-NFRLength.bed
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

TF=ATF7
JASPAR=MA0834-1
ENCFF=ENCFF868QLL
ODIR=Intersect/$TF

# Intersect with NFR RefPT for NFR overlap or not (a.k.a. "promoter proximal" or not)
bedtools intersect -wa -u -a $ODIR/BoundMotifs.bed -b $NFR_BEDFILE > $MOTIF/$TF\_Bound-Promoter.bed
bedtools intersect -wa -v -a $ODIR/BoundMotifs.bed -b $NFR_BEDFILE > $MOTIF/$TF\_Bound-NonPromoter.bed
# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_Bound-Promoter.bed -o $MOTIF/1000bp/$TF\_Bound-Promoter_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_Bound-NonPromoter.bed -o $MOTIF/1000bp/$TF\_Bound-NonPromoter_1000bp.bed
# Print line count stats
wc -l $ODIR/BoundMotifs.bed
wc -l $MOTIF/$TF\_Bound-Promoter.bed
wc -l $MOTIF/$TF\_Bound-NonPromoter.bed
