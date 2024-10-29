#!/bin/bash

# Make PE insert size histogram of BNase-seq, MNase-seq ([21 U],[304 U]), and DNase-seq
# (1c) BNase-seq
# (1e) MNase-seq
# (1f) DNase-seq
# see 04_Figures/Fig_1c.sh
# see 04_Figures/Fig_1e_f.sh

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/Z_Figures
###

# Dependencies
# - java
# - pandas
# - python
# - seaborn

set -exo
module load anaconda3
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
HISTOGRAM=$WRK/../bin/make_fragment_histograms.py

[ -d F1/c ] || mkdir -p F1/c

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_merge_hg38.bam -o F1/c/BNase-seq_50U-10min_merge_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/DNase-seq_ENCFF518XTC_rep1_hg38.bam -o F1/c/DNase-seq_ENCFF518XTC_rep1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/MNase-seq_21U_rep1_hg38.bam -o F1/c/MNase-seq_21U_rep1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/MNase-seq_304U_rep1_hg38.bam -o F1/c/MNase-seq_304U_rep1_hg38

# Generate insert size frequency histograms
python $HISTOGRAM --ymax 16000000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160 -i F1/c/BNase-seq_50U-10min_merge_hg38_InsertHistogram.out -o F1/c/BNase-seq_50U-10min_merge_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax  1200000 --reflines 70 80 90 100 110 120 130 140 -i F1/c/DNase-seq_ENCFF518XTC_rep1_hg38_InsertHistogram.out -o F1/c/DNase-seq_ENCFF518XTC_rep1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax  1400000 --reflines 106 113 128 134 149 -i F1/c/MNase-seq_21U_rep1_hg38_InsertHistogram.out -o F1/c/MNase-seq_21U_rep1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax  1000000 --reflines 106 127 145 -i F1/c/MNase-seq_304U_rep1_hg38_InsertHistogram.out -o F1/c/MNase-seq_304U_rep1_hg38_InsertHistogram_HIST.svg