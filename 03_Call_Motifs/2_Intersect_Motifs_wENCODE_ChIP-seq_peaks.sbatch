#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/2_Intersect_Motifs_wENCODE_ChIP.log.out-%a
#SBATCH -e logs/2_Intersect_Motifs_wENCODE_ChIP.log.err-%a
#SBATCH --array 1-11

# FIMO the reference genome for each motif in the PWM directory

### CHANGE ME
METADATA=TF_JASPAR_ENCODE_config.txt
WRK=/path/to/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/03_Call_Motifs
###

# Dependencies
# - bedtools
# - java
# - MEME suite (FIMO)

set -exo
module load anaconda3
module load bedtools
source activate bx

# Inputs and outputs
GENOME=$WRK/../data/hg19_files/hg19.fa
BLACKLIST=$WRK/../data/hg19_files/hg19_exclude.bed
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
# Intersect peaks with motifs
bedtools intersect -wo -a $ODIR/ENCODE_shuffled_1000bp.bed -b $ODIR/JASPAR_shuffled_20bp.bed > $ODIR/Intersect.bed
# Count sites and divide by two (# lower half)
NSITES=`wc -l $ODIR/Intersect.bed | awk '{print int($1/2)}'`
echo NSITES=$NSITES
# cat $TARGET_INTERSECT | wc -l | awk '{printf \"%.f\\n\", \$1 * 0.5}' > $NUMBER
# Sort by ChIP enrichment value ("signal value") column 7 (from ENCODE peak file), then make bedfile from columns 11-16 (from motif bedfile).
sort -rnk7,7 $ODIR/Intersect.bed | tail -n $NSITES | \
	awk '{FS="\t"}{OFS="\t"}{print $11,$12,$13,$14,$15,$16}' > Intersect/$TF/LowerBound.bed

[ -d $MOTIF ] || mkdir $MOTIF
[ -d $MOTIF/1000bp ] || mkdir $MOTIF/1000bp

mv Intersect/$TF/LowerBound.bed $MOTIF/$TF\_LowerBound.bed
# Expand 1000bp (Motif-centered lower-bound intersect)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MOTIF/$TF\_LowerBound.bed -o $MOTIF/1000bp/$TF\_LowerBound_1000bp.bed
