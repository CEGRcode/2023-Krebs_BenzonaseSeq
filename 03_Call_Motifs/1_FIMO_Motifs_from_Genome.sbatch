#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 01:00:00
#SBATCH -A open
#SBATCH -o logs/1_FIMO_Motifs_from_Genome.log.out-%a
#SBATCH -e logs/1_FIMO_Motifs_from_Genome.log.err-%a
#SBATCH --array 1-89

# FIMO the reference genome for each motif in the PWM directory

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/03_Call_Motifs
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/03_Call_Motifs
METADATA=TF_JASPAR_ENCODE_config.txt
###

# Dependencies
# - bedtools
# - java
# - MEME suite (FIMO)

set -exo
module load bedtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BLACKLIST=$WRK/../data/hg38_files/ENCFF356LFX_hg38_exclude.bed

# Script shortcuts
DEDUP=$WRK/../bin/dedup_coord_by_ID.py
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Determine TF name, JASPAR motif name, and ENCODE peak file name
TF=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
JASPAR=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
ENCFF=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Construct input filenames
JDIR=$WRK/../data/JASPAR
PWMFILE=$JDIR/$TF\_$JASPAR.meme
FIMO=FIMO/$TF\_$JASPAR

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d FIMO ] || mkdir FIMO
[ -d $FIMO ] || mkdir $FIMO

echo "($SLURM_ARRAY_TASK_ID) $TF"

# Run FIMO to scan genome for motif occurrences
# fimo --verbosity 1 --thresh 1.0E-4  --max-strand --bfile --motif-- --oc FIMO_original $PWMFILE $GENOME
fimo --verbosity 1 --thresh 1.0E-4 --max-strand         --oc $FIMO $PWMFILE $GENOME

# Convert GFF to BED format
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed <(sed '1d' $FIMO/fimo.gff) -o $FIMO/fimo_unformatted_unfiltered.bed

# Rename motif instances by coord info and sort by id, then FIMO motif score
awk '{OFS="\t"}{FS="\t"}{print $1,$2,$3,"M-"$1"_"$2"_"$3,$5,$6}' $FIMO/fimo_unformatted_unfiltered.bed \
	| sort -k4 -rnk5 > $FIMO/fimo_unfiltered.bed

# Handle palindromic motifs with dedup keeping higher FIMO score
python $DEDUP -i $FIMO/fimo_unfiltered.bed -o $FIMO/fimo.bed

# Expand 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $FIMO/fimo.bed -o $FIMO/fimo_1000bp.bed

# Filter blacklist (exclude coords expanded off chromosome)
bedtools intersect -v -a <(awk '{if($2>0) print}' $FIMO/fimo_1000bp.bed) -b $BLACKLIST > $FIMO/filtered.bed
