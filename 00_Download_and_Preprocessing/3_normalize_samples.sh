#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=14gb
#PBS -l walltime=00:10:00
#PBS -A open
#PBS -o logs/get_scaling_factors.log.out
#PBS -e logs/get_scaling_factors.log.err
#PBS -t 1-NSAMPLES

# Calculate normalization factor for every *.bam file in a directory

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq
BAMDIR=$WRK/data/BAM
BEDDIR=$WRK/data/BED
FDIR=$WRK/data/NormalizationFactors
BLACKLIST=$WRK/data/ChExMix_Peak_Filter_List_190612.bed
CONTROL=$WRK/data/BAM/masterNoTag_20180928.bam
#PBS_ARRAYID=1
###

module load gcc
module load samtools
module load anaconda3
#source activate my-env

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14.jar
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.14-$PBS_ARRAYID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
# (none)

# Set up output directories
cd $WRK
[ -d logs ] || mkdir logs
[ -d $FDIR ] || mkdir $FDIR

# Determine BAM file for the current job array index
BAMFILE=`ls $BAMDIR/*.bam | head -n $PBS_ARRAYID | tail -1`;
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Use different normalization method depending on target/assay
# This script ONLY covers NCIS (for TF ChIPs), if your data includes NUCLEOSOME data, USE NFR normalization!
echo "Calculate classic TF NCIS normalization factors w/ blacklist"
java -jar $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE -f $BLACKLIST --ncis -c $CONTROL -w 500 -o $WRK/data/NormalizationFactors/$BAM\_NCISb

rm $SCRIPTMANAGER
