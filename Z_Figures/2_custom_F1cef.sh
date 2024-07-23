#!/bin/bash

# Make PE insert size histogram of...
# (1c) BNase-seq
# (1e) MNase-seq
# (1f) DNase-seq
# see 04_Figures/Fig_1c.sh
# see 04_Figures/Fig_1e_f.sh

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/Z_Figures
###

# Dependencies
# - java
# - pandas
# - python
# - seaborn

set -exo
module load anaconda3
source activate bx

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
HISTOGRAM=$WRK/../bin/make_fragment_histograms-MOD.py
HISTOGRAM2=$WRK/../bin/make_fragment_histograms-MOD2.py

[ -d F1/c ] || mkdir -p F1/c
[ -d F1/e ] || mkdir -p F1/e
[ -d F1/f ] || mkdir -p F1/f

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U_merge_hg19.bam -o F1/c/BNase-seq_50U_merge_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/MNase-seq_21U_rep1_hg19.bam -o F1/e/MNase-seq_21U_rep1_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/MNase-seq_304U_rep1_hg19.bam -o F1/e/MNase-seq_304U_rep1_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/DNase-seq_-_rep1_hg19.bam -o F1/f/DNase-seq_-_rep1_hg19

# Generate insert size frequency histograms
python $HISTOGRAM -i F1/c/BNase-seq_50U_merge_hg19_InsertHistogram.out -o F1/c/BNase-seq_50U_merge_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM2 -i F1/c/BNase-seq_50U_merge_hg19_InsertHistogram.out -o F1/c/BNase-seq_50U_merge_hg19_InsertHistogram_HIST2.svg
python $HISTOGRAM -i F1/e/MNase-seq_21U_rep1_hg19_InsertHistogram.out -o F1/e/MNase-seq_21U_rep1_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM -i F1/e/MNase-seq_304U_rep1_hg19_InsertHistogram.out -o F1/e/MNase-seq_304U_rep1_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM -i F1/f/DNase-seq_-_rep1_hg19_InsertHistogram.out -o F1/f/DNase-seq_-_rep1_hg19_InsertHistogram_HIST.svg
