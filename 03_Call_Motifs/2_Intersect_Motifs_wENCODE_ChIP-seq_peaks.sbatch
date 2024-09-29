#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/2_Intersect_Motifs_wENCODE_ChIP.log.out-%a
#SBATCH -e logs/2_Intersect_Motifs_wENCODE_ChIP.log.err-%a
#SBATCH --array 1-11

# FIMO the reference genome for each motif in the PWM directory

### CHANGE ME
METADATA=TF_JASPAR_ENCODE_config.txt
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
GENOME=$WRK/../data/hg38_files/hg38.fa
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed

# Script shortcuts
DEDUP=$WRK/../bin/dedup_coord_by_ID.py
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Determine TF name, JASPAR motif name, and ENCODE peak file name
TF=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
FJASPAR=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}' | sed -e 's/\./-/g'`
ENCFF=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Construct input filenames
JDIR=$WRK/../data/JASPAR
MOTIF=$WRK/../data/RefPT-Motif
PWMFILE=$JDIR/$TF\_$FJASPAR.meme

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d Intersect ] || mkdir Intersect
[ -d Intersect/$TF ] || mkdir Intersect/$TF

echo "($SLURM_ARRAY_TASK_ID) $TF"

ODIR=Intersect/$TF

# Shuffle motif instances and ENCODE peaks (avoid ordering effects)
shuf FIMO/$TF/filtered.bed > $ODIR/JASPAR_shuffled.bed
gzip -dc narrowPeaks/$TF\_$ENCFF.bed.gz | shuf > $ODIR/ENCODE_shuffled.bed
# Expand JASPAR by 20bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 20   $ODIR/JASPAR_shuffled.bed -o $ODIR/JASPAR_shuffled_20bp.bed
# Expand ENCODE by 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $ODIR/ENCODE_shuffled.bed -o $ODIR/ENCODE_shuffled_1000bp.bed
# Intersect peaks with motifs - filter to keep overlap - move ENCODE ChIP value ("signal value") to score col - sort by ID, then score
bedtools intersect -loj -a Intersect/$TF/JASPAR_shuffled_20bp.bed -b Intersect/$TF/ENCODE_shuffled_1000bp.bed \
	| awk '{OFS="\t"}{FS="\t"}{if($8>0) print $1,$2,$3,$4,$13,$6}' | sort -rnk4,5 \
	> $ODIR/BoundMotifs-wDup.bed
# Deduplicate motifs by keeping first instance (larger ENCODE score based on previous command sort)
python $DEDUP -i $ODIR/BoundMotifs-wDup.bed -o $ODIR/BoundMotifs.bed
# Count sites and divide by two (# lower half)
NSITES=`wc -l $ODIR/BoundMotifs.bed | awk '{print int($1/2)}'`
echo NSITES=$NSITES
# Sort by ChIP enrichment and take bottom half
sort -rnk5,5 $ODIR/BoundMotifs.bed | tail -n $NSITES > Intersect/$TF/LowerBound.bed

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

mv Intersect/$TF/LowerBound.bed $MOTIF/$TF\_Bound-LowerHalf.bed
# Expand 1000bp (Motif-centered lower-bound intersect)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_Bound-LowerHalf.bed -o $MOTIF/1000bp/$TF\_Bound-LowerHalf_1000bp.bed
