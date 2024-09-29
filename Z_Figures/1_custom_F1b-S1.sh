#!/bin/bash

# (F1b and S1) Use peak-align to "pileup" CpG island annotations on TSS-centered RefPT
# (F1b) Tag Pileup deep ENCODE MNase-seq signal on TSS-centered RefPT

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/Z_Figures
THREADS=4
###

# Dependencies
# - java

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Set up output directories
[ -d F1/b ] || mkdir -p F1/b
[ -d S1 ] || mkdir S1

REFPT=$WRK/../data/RefPT-Krebs/2000bp/TSS_GROUP-Expressed_SORT-CpG_2000bp.bed
PEAK=$WRK/../data/RefPT-Other/CpG_Islands.bed
BAMFILE=$WRK/../data/BAM/MNase-seq_ENCODE_merge_hg19.bam

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.14.jar
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl

BED=`basename $REFPT ".bed"`
BAM=`basename $BAMFILE ".bam"`

# ===============================================================================================================================

echo "Peak-align CpG islands around TSS"
BASE=CpG-Islands_$BED

# Pileup CpG islands
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref $PEAK $REFPT -o F1/b/$BASE.cdt

# Two-color heatmap
java -jar $SCRIPTMANAGER figure-generation heatmap -a 1 --blue F1/b/$BASE.cdt -o F1/b/$BASE\_treeview.png

# Count sites
NSITES=`wc -l $REFPT | awk '{print $1-1}'`

# Label heatmap
java -jar $SCRIPTMANAGER figure-generation label-heatmap F1/b/$BASE\_treeview.png \
	-l "-1" -m "0" -r "+1" -w 2 -f 18 \
	-x "Distance from TSS (kb)" -y "${NSITES} CoPRO determined TSSs sorted by CpG island length" \
	-o F1/b/$BASE\_treeview.svg

# Copy same file into supplement
cp F1/b/$BASE\_treeview.svg S1/$BASE\_treeview.svg

# ===============================================================================================================================

echo "Pileup MNase-seq across TSS RefPT w/ 80bp shift"
BASE=$BAM\_$BED

# Pileup SE MNase data (shift 80bp)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $REFPT $BAMFILE \
	-5 -1 --shift 80 --combined --cpu $THREADS -z \
	-o F1/b/$BASE\_combined.out -M F1/b/$BASE

# Two-color heatmap
java -jar $SCRIPTMANAGER figure-generation heatmap --black -p 0.95 F1/b/$BASE\_combined.cdt.gz -o F1/b/$BASE\_treeview.png

# Label heatmap
java -jar $SCRIPTMANAGER figure-generation label-heatmap F1/b/$BASE\_treeview.png \
	-l "-1" -m "0" -r "+1" -w 2 -f 18 \
	-x "Distance from TSS (kb)" -y "${NSITES} CoPRO determined TSSs sorted by CpG island length" \
	-o F1/b/$BASE\_treeview.svg
