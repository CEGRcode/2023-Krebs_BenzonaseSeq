#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/2_Intersect_Motifs_wENCODE_ChIP.log.out-%a
#SBATCH -e logs/2_Intersect_Motifs_wENCODE_ChIP.log.err-%a
#SBATCH --array 1-11

# Re-intersect all ATF7 motif with ChIP-seq peaks and sort into
# promoter-proximal and not promoter-proximal (NFR overlap or not)

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/03_Call_Motifs
###

# Dependencies
# - bedtools
# - javas

set -exo
module load anaconda3
module load bedtools
source activate bx

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif
NFR_BEDFILE=$WRK/../data/RefPT-Krebs/NFR_K562.bed
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar

TF=ATF7
JASPAR=MA0834-1
ENCFF=ENCFF868QLL
ODIR=Intersect/$TF

# Expand ENCODE by 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $ODIR/ENCODE_shuffled.bed -o $ODIR/ENCODE_shuffled_1000bp.bed
# Intersect JASPAR motif instances with ENCODE ChIP-seq peaks and write the motif-centered hits
bedtools intersect -wb -a $ODIR/ENCODE_shuffled_1000bp.bed -b $ODIR/JASPAR_shuffled_20bp.bed > $ODIR/Intersect_Motif-centered.bed
# Intersect with NFR RefPT for NFR overlap or not (a.k.a. "promoter proximal" or not)
bedtools intersect -u -a $ODIR/Intersect_Motif-centered.bed -b $NFR_BEDFILE > $MOTIF/$TF\_PromoterProximal.bed
bedtools intersect -v -a $ODIR/Intersect_Motif-centered.bed -b $NFR_BEDFILE > $MOTIF/$TF\_NotPromoterProximal.bed
# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_PromoterProximal.bed -o $MOTIF/1000bp/$TF\_PromoterProximal_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_NotPromoterProximal.bed -o $MOTIF/1000bp/$TF\_NotPromoterProximal_1000bp.bed
# Print line count stats
wc -l $ODIR/Intersect_Motif-centered.bed
wc -l $MOTIF/$TF\_PromoterProximal.bed
wc -l $MOTIF/$TF\_NotPromoterProximal.bed
