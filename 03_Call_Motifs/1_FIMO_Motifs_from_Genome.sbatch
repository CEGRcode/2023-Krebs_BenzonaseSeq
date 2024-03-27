#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 20:00:00
#SBATCH -A open
#SBATCH -o logs/1_FIMO_Motifs_from_Genome.log.out-%a
#SBATCH -e logs/1_FIMO_Motifs_from_Genome.log.err-%a
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
PWMFILE=$JDIR/$TF\_$FJASPAR.meme

# Create output directories if they don't exist
[ -d logs ] || mkdir logs
[ -d FIMO_default ] || mkdir FIMO_default
[ -d FIMO ] || mkdir FIMO
[ -d FIMO/$TF ] || mkdir FIMO/$TF

echo "($SLURM_ARRAY_TASK_ID) $TF"

# Run FIMO to scan genome for motif occurrences
#fimo --oc FIMO_default/$TF           --verbosity 1 --thresh 1.0E-4 $PWMFILE $GENOME
#fimo --oc FIMO/$TF --bfile --motif-- --verbosity 1 --thresh 1.0E-4 $PWMFILE $GENOME

# Convert GFF to BED format
java -jar $SCRIPTMANAGER coordinate-manipulation gff-to-bed <(sed '1d' FIMO/$TF/fimo.gff) -o FIMO/$TF/fimo.bed

# Expand 1000 bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 FIMO/$TF/fimo.bed -o FIMO/$TF/fimo_1000bp.bed

# Filter blacklist (exclude coords expanded off chromosome)
bedtools intersect -v -a <(awk '{if($2>0) print}' FIMO/$TF/fimo_1000bp.bed) -b $BLACKLIST > FIMO/$TF/filtered.bed
