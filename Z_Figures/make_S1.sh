#!/bin/bash

# Make PE insert size histogram of BNase-seq 2 replicates at each of three digestion levels (50U-3min, 50U-10min, 50U-30min)

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

[ -d S1 ] || mkdir -p S1

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_1_hg38.bam -o S1/BNase-seq_50U-3min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_2_hg38.bam -o S1/BNase-seq_50U-3min_2_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_1_hg38.bam -o S1/BNase-seq_50U-10min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_2_hg38.bam -o S1/BNase-seq_50U-10min_2_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_1_hg38.bam -o S1/BNase-seq_50U-30min_1_hg38
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_2_hg38.bam -o S1/BNase-seq_50U-30min_2_hg38

# Generate insert size frequency histograms
python $HISTOGRAM --ymax 100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/BNase-seq_50U-3min_1_hg38_InsertHistogram.out -o S1/BNase-seq_50U-3min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/BNase-seq_50U-3min_2_hg38_InsertHistogram.out -o S1/BNase-seq_50U-3min_2_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 400000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/BNase-seq_50U-10min_1_hg38_InsertHistogram.out -o S1/BNase-seq_50U-10min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 5100000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/BNase-seq_50U-10min_2_hg38_InsertHistogram.out -o S1/BNase-seq_50U-10min_2_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 140000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/BNase-seq_50U-30min_1_hg38_InsertHistogram.out -o S1/BNase-seq_50U-30min_1_hg38_InsertHistogram_HIST.svg
python $HISTOGRAM --ymax 140000 --reflines 30 40 50 60 70 80 90 100 110 120 130 140 150 160  -i S1/BNase-seq_50U-30min_2_hg38_InsertHistogram.out -o S1/BNase-seq_50U-30min_2_hg38_InsertHistogram_HIST.svg
