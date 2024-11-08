#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/2_Intersect_Motifs_wENCODE_ChIP-seq_peaks.log.out-%a
#SBATCH -e logs/2_Intersect_Motifs_wENCODE_ChIP-seq_peaks.log.err-%a
#SBATCH --array 1-89

# FIMO the reference genome for each motif in the PWM directory
# see 241017_TF_scripts/01_TF_pipeline/TF_pipeline_v13_241011.sh

# data/RefPT-JASPAR
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile1.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile2.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile3.bed
#   |--<TFNAME>_<JASPARID>_SORT-TFnucRatio_GROUP-Quartile4.bed

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/03_Call_Motifs
METADATA=TF_JASPAR_ENCODE_config.txt
###

# Dependencies
# - bedtools
# - java

set -exo
module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BAMFILE=$WRK/../data/BAM/BNase-seq_50U-10min_merge_hg38.bam
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed
MOTIF=$WRK/../data/RefPT-JASPAR

# Script shortcuts
DEDUP=$WRK/../bin/dedup_coord_by_ID.py
RATIO=$WRK/../bin/calculate_BED_ScoreRatio.pl
UPDATES=$WRK/../bin/update_BED_score_with_TAB_score.pl
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Determine TF name, JASPAR motif name, and ENCODE peak file name
TF=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
JASPAR=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
ENCFF=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Construct more input output files and directories
JDIR=$WRK/../data/JASPAR
PWMFILE=$JDIR/$TF\_$JASPAR.meme
FIMO=FIMO/$TF\_$JASPAR
ODIR=Intersect/$TF\_$JASPAR

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d Intersect ] || mkdir Intersect
[ -d $ODIR ] || mkdir $ODIR

echo "($SLURM_ARRAY_TASK_ID) $TF"

## =====Associate motifs with ENCODE peaks=====

# Parse BED 6 + 3 narroPeak files into summits (1bp)
gzip -dc narrowPeaks/$TF\_$ENCFF.bed.gz \
	| awk '{OFS="\t"}{FS="\t"}{print $1,$2+$10,$2+$10+1,"P-"$1"_"$2"_"$3,$5,$6}' \
	> $ODIR/ENCODE_summit.bed
# Read more on narrowPeak format here: https://genome.ucsc.edu/FAQ/FAQformat.html#format12

# Sort motif instances and ENCODE peaks (for closest command later)
bedtools sort -i $FIMO/filtered.bed > $ODIR/JASPAR.bed
bedtools sort -i $ODIR/ENCODE_summit.bed > $ODIR/ENCODE.bed

# Get closest motif for each peak up to 1000bp away
bedtools closest -d -a $ODIR/ENCODE.bed -b $ODIR/JASPAR.bed \
	| awk '{OFS="\t"}{FS="\t"}{if ($13<1000) print}' \
	> $ODIR/ENCODE_matchedToJASPAR.tsv

# Select closest peak-motif match when multiple motifs match a peak
sort -nk13,13 $ODIR/ENCODE_matchedToJASPAR.tsv \
	| sort -uk4,4 > $ODIR/UniquePeakMatch.tsv

# Identify unique peak-motif match when multiple peaks match a motif (pick greater ENCODE binding score, then best FIMO score)
sort -nk7,7 -nk11,11 $ODIR/UniquePeakMatch.tsv \
	| sort -uk10,10 \
	> $ODIR/UniqueMatch.tsv

## =====Sort by TFsignal : Nucleosome occupancy ratios=====

# Reformat as motif BED with ENCODE ChIP enrichment for score
awk 'BEGIN{OFS="\t";FS="\t"}{print $7,$8,$9,$10,$5,$12}' $ODIR/UniqueMatch.tsv > $ODIR/BoundMotifs_SCORE-TFsignal.bed

# Expand 164bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 164 $ODIR/BoundMotifs_SCORE-TFsignal.bed -o $ODIR/BoundMotifs_SCORE-TFsignal_164bp.bed

# Tag Pileup to calculate Nucleosome occupancy score
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m --combined $ODIR/BoundMotifs_SCORE-TFsignal_164bp.bed $BAMFILE -M $ODIR/BNase_BoundMotif

# Sum scores to a per-ID value
java -jar $SCRIPTMANAGER read-analysis aggregate-data $ODIR/BNase_BoundMotif_combined.cdt -o $ODIR/BNase_BoundMotif_NucleosomeOccupancy.tab

# Add nucleosome occupancy score to BED file
perl $UPDATES $ODIR/BoundMotifs_SCORE-TFsignal.bed $ODIR/BNase_BoundMotif_NucleosomeOccupancy.tab $ODIR/BoundMotifs_SCORE-NucleosomeOccuppancy.bed

# Calculate ratio
perl $RATIO $ODIR/BoundMotifs_SCORE-TFsignal.bed $ODIR/BoundMotifs_SCORE-NucleosomeOccuppancy.bed $ODIR/BoundMotifs_SCORE-TFnucRatio.bed

# Sort by ratio
sort -nk5,5 $ODIR/BoundMotifs_SCORE-TFnucRatio.bed > $ODIR/BoundMotifs_SORT-TFnucRatio.bed

## =====Split into Quartiles=====

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

# Count sites and divide by four (quartile size)
NSITES=`wc -l $ODIR/BoundMotifs_SORT-TFnucRatio.bed | awk '{print int($1/4)}'`
TWON=`wc -l $ODIR/BoundMotifs_SORT-TFnucRatio.bed | awk '{print int($1*2/4)}'`
THREEN=`wc -l $ODIR/BoundMotifs_SORT-TFnucRatio.bed | awk '{print int($1*3/4)}'`
echo NSITES=$NSITES

# Slice quartiles from counts
head -n $NSITES $ODIR/BoundMotifs_SORT-TFnucRatio.bed                > $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile1.bed
head -n $TWON   $ODIR/BoundMotifs_SORT-TFnucRatio.bed | tail -n $NSITES > $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile2.bed
head -n $THREEN $ODIR/BoundMotifs_SORT-TFnucRatio.bed | tail -n $NSITES > $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile3.bed
tail -n $NSITES $ODIR/BoundMotifs_SORT-TFnucRatio.bed                > $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile4.bed

# Expand 1000bp (Motif-centered lower-bound intersect)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile1.bed -o $MOTIF/1000bp/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile1_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile2.bed -o $MOTIF/1000bp/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile2_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile3.bed -o $MOTIF/1000bp/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile3_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile4.bed -o $MOTIF/1000bp/$TF\_$JASPAR\_SORT-TFnucRatio_GROUP-Quartile4_1000bp.bed
