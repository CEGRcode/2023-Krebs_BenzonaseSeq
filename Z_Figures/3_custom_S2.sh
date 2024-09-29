#!/bin/bash

# (S2) Make PE insert size histogram of BNase-seq titrations.
# see 04_Figures/ED_Fig_2.sh

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
source activate bx

# Inputs and outputs
BAMDIR=$WRK/../data/BAM

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
HISTOGRAM=$WRK/../bin/make_fragment_histograms-MOD.py
HISTOGRAM2=$WRK/../bin/make_fragment_histograms-MOD2.py

[ -d S2 ] || mkdir S2

# "28385_Input_-_K562_raw_master.bam" "28386_Input_-_K562_raw_master.bam" "28389_Input_-_K562_raw_master.bam" "28390_Input_-_K562_raw_master.bam" "28393_Input_-_K562_raw_master.bam" "28394_Input_-_K562_raw_master.bam";

# Run Paired-end Statistics
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_rep1_hg19.bam -o S2/BNase-seq_50U-3min_rep1_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-3min_rep2_hg19.bam -o S2/BNase-seq_50U-3min_rep2_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_rep1_hg19.bam -o S2/BNase-seq_50U-10min_rep1_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-10min_rep2_hg19.bam -o S2/BNase-seq_50U-10min_rep2_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_rep1_hg19.bam -o S2/BNase-seq_50U-30min_rep1_hg19
java -jar $SCRIPTMANAGER bam-statistics pe-stat --min 0 --max 340 $BAMDIR/BNase-seq_50U-30min_rep2_hg19.bam -o S2/BNase-seq_50U-30min_rep2_hg19

# Generate insert size frequency histograms
python $HISTOGRAM -i S2/BNase-seq_50U-3min_rep1_hg19_InsertHistogram.out -o S2/BNase-seq_50U-3min_rep1_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM -i S2/BNase-seq_50U-3min_rep2_hg19_InsertHistogram.out -o S2/BNase-seq_50U-3min_rep2_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM -i S2/BNase-seq_50U-10min_rep1_hg19_InsertHistogram.out -o S2/BNase-seq_50U-10min_rep1_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM -i S2/BNase-seq_50U-10min_rep2_hg19_InsertHistogram.out -o S2/BNase-seq_50U-10min_rep2_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM -i S2/BNase-seq_50U-30min_rep1_hg19_InsertHistogram.out -o S2/BNase-seq_50U-30min_rep1_hg19_InsertHistogram_HIST.svg
python $HISTOGRAM -i S2/BNase-seq_50U-30min_rep2_hg19_InsertHistogram.out -o S2/BNase-seq_50U-30min_rep2_hg19_InsertHistogram_HIST.svg
