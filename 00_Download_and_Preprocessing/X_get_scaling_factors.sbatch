#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:10:00
#SBATCH -A open
#SBATCH -o logs/get_scaling_factors.log.out-%a
#SBATCH -e logs/get_scaling_factors.log.err-%a
#SBATCH --array 1-NSAMPLES

# NSAMPLES should be replaced with `ls ../data/BAM/*.bam | wc -l`

# SLURM_ARRAY_TASK_ID=2

# Calculate normalization factor for every *.bam file in a directory

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing/
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
###

# Dependencies
# - java
# - python
# - pandas
# - seaborn

set -exo
module load samtools
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BAMDIR=$WRK/../data/BAM
FDIR=$WRK/../data/BAM/NormalizationFactors
BLACKLIST=$WRK/../data/hg38_files/hg38_exclude.bed
# BLACKLIST=$WRK/../data/hg19_files/hg19_exclude.bed

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
# (none)

# Set up output directories
cd $WRK
[ -d logs ] || mkdir logs
[ -d $FDIR ] || mkdir $FDIR

# Determine BAM file for the current job array index
BAMFILE=`ls $BAMDIR/*.bam | head -n $SLURM_ARRAY_TASK_ID | tail -1`;
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Use different normalization method depending on target/assay
echo "Calculate Total Tag normalization factors w/ blacklist"
java -jar $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE -f $BLACKLIST --total-tag -o $FDIR/$BAM\_TotalTag

rm $SCRIPTMANAGER
