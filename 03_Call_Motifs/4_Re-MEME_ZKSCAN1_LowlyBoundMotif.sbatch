#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 01:00:00
#SBATCH -A open
#SBATCH -o logs/4_Re-MEME_ZKSCAN1_LowlyBoundMotif.log.out
#SBATCH -e logs/4_Re-MEME_ZKSCAN1_LowlyBoundMotif.log.err

# Re-MEME off lowly bound ZKSCAN1 sites to get an extended motif and tally
# the number of bound motifs from the full set and the lower-half that contaiin
# the extended motif.

# Note: As of some 20240404 checks, it seems the extended motif does not
# specifically reflect the "Lowly Bound" sites  as it occurs across most bound
# instances. May remove this later.

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/03_Call_Motifs
###

# Dependencies
# - bedtools
# - java
# - MEME suite v5.0.5 (MEME,FIMO,meme-get-motif)

set -exo
module load anaconda3
module load bedtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
MOTIF=$WRK/../data/RefPT-Motif
NFR_BEDFILE=$WRK/../data/RefPT-Krebs/NFR_K562.bed
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
GENOME=$WRK/../data/hg38_files/hg38.fa
OUTPUT=ZKSCAN1_ExtendedMotif

TF=ZKSCAN1
JASPAR=MA1585-1
ENCFF=ENCFF163VUK

[ -d $OUTPUT ] || mkdir $OUTPUT
[ -d $OUTPUT/FIMO-Full ] || mkdir $OUTPUT/FIMO-Full
[ -d $OUTPUT/FIMO-Lower ] || mkdir $OUTPUT/FIMO-Lower
[ -d $OUTPUT/MEME ] || mkdir $OUTPUT/MEME

# Extract FASTA (Lowly bound)
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $MOTIF/1000bp/$TF\_Bound-LowerHalf_1000bp.bed -o $OUTPUT/$TF\_Bound-LowerHalf_1000bp.fa
# Re-MEME for lowly bound motif
meme -dna -oc $OUTPUT/MEME -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 100 \
	-objfun classic -revcomp -markov_order 0 $OUTPUT/$TF\_Bound-LowerHalf_1000bp.fa
# Extract Motif 1 from MEME.txt file (5.0.5)
MEMEID=`grep "MOTIF [A-Z]\+ MEME-1" $OUTPUT/MEME/meme.txt | awk '{print $2}'`
meme-get-motif -id $MEMEID $OUTPUT/MEME/meme.txt > $OUTPUT/$TF\_Bound-LowerHalf_Re-MEME-M1.meme.txt

# Expand 1000bp (Full Bound set)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 Intersect/$TF/BoundMotifs.bed -o Intersect/$TF/BoundMotifs_1000bp.bed
# Extract FASTA (Full Bound set)
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME Intersect/$TF/BoundMotifs_1000bp.bed -o Intersect/$TF/BoundMotifs_1000bp.fa

# FIMO lowly bound for instances of the extended motif
fimo --oc $OUTPUT/FIMO-Lower --verbosity 1 --bfile --nrdb-- --thresh 1.0E-4 --norc $OUTPUT/$TF\_Bound-LowerHalf_Re-MEME-M1.meme.txt $OUTPUT/$TF\_Bound-LowerHalf_1000bp.fa
# FIMO all bound sites for instances of the extended motif
fimo --oc $OUTPUT/FIMO-Full  --verbosity 1 --bfile --nrdb-- --thresh 1.0E-4 --norc $OUTPUT/$TF\_Bound-LowerHalf_Re-MEME-M1.meme.txt Intersect/$TF/BoundMotifs_1000bp.fa

NSITES_LOWER=`wc -l $OUTPUT/$TF\_Bound-LowerHalf_1000bp.fa | awk '{print $1/2}'`
NSITES_FULL=`wc -l Intersect/$TF/BoundMotifs_1000bp.fa | awk '{print $1/2}'`
NSITES_LOWER_FIMO=`cut -f1 $OUTPUT/FIMO-Lower/fimo.gff | sort | uniq | wc -l | awk '{print $1-1}'`
NSITES_FULL_FIMO=`cut -f1 $OUTPUT/FIMO-Full/fimo.gff | sort | uniq | wc -l | awk '{print $1-1}'`

echo "LowerHalf of Bound: $NSITES_LOWER_FIMO FIMO'd of $NSITES_LOWER"
echo "Full set of Bound: $NSITES_FULL_FIMO FIMO'd of $NSITES_FULL"
