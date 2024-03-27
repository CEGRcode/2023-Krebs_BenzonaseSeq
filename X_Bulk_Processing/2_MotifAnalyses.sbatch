#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/2_MotifAnalyses.log.out-%a
#SBATCH -e logs/2_MotifAnalyses.log.err-%a
#SBATCH --array 1-11

# Loop through each BEDFILE in data/RefPT-Motif/1000bp and perform motif analyses:
# - BI pileups

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/X_Bulk_Processing
METADATA=BI_Pileups.txt
THREADS=4
###

# Dependencies
# - java

set -exo
module load anaconda
conda activate bx

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.14.jar
RESIZE=../bin/resize_png.py
DINUCLEOTIDE=../bin/
COMPOSITE=../bin/sum_Col_CDT.pl

# Inputs and outputs
GENOME=$WRK/../data/hg19_files/hg19.fa
MOTIF=$WRK/../data/RefPT-Motif
OUTDIR=Library

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

BAMFILE=$WRK/../data/BAM/BNase-seq_50U_merge_hg19.bam
BAM=`basename $BAMFILE ".bam"`

# Determine BED file for the current job array index
BEDFILE=`ls $MOTIF/1000bp/*.bed | head -n $SLURM_ARRAY_TASK_ID | tail -1`
BED=`basename $BEDFILE ".bed"`

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1}'`

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} "x" ${BAMFILE} "

DIR=$OUTDIR/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

# ===============================================================================================================================

echo "Run Motif endo cuts BI pileup"
BASE=$BAM\_$BED\_5both

# Pileup (endo)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu $THREADS -5 -a -o $DIR/Composites/$BASE.out
# No scaling
