#!/bin/bash

# Make PE insert size histogram of BNase-seq, MNase-seq ([21 U],[304 U]), and DNase-seq
# (1c) BNase-seq
# (1e) MNase-seq
# (1f) DNase-seq
# see 04_Figures/Fig_1c.sh
# see 04_Figures/Fig_1e_f.sh
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
# - pandas
# - python
# - seaborn

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Inputs and outputs
KREBS=$WRK/../data/RefPT-Krebs/
REFPT=$KREBS/2000bp/TSS_GROUP-Expressed_SORT-CpG_2000bp.bed
LIBRARY=$WRK/../X_Bulk_Processing/Library

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl
VIOLIN=$WRK/../bin/make_violin_plot.py

# Set up output directories
[ -d F7/a ] || mkdir -p F7/a
BED=`basename $REFPT ".bed"`

# ===============================================================================================================================

echo "Copy heatmaps over from X_Bulk_Processing/Library"

cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/BNase-seq_50U-10min_merge_hg38_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg F7/a/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/DNase-seq_ENCFF518XTC_rep1_hg38_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg F7/a/


# ===============================================================================================================================

CPGBED=$WRK/../data/RefPT-Other/CpGIslands.bed

echo "Peak-align CpG islands around TSS"
BASE=CpG-Islands_$BED

# Pileup CpG islands
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref $CPGBED $REFPT -o F7/a/$BASE

# Two-color heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -a 1 --blue F7/a/$BASE\_combined.cdt -o F7/a/$BASE\_treeview.png

# Count sites
NSITES=`wc -l $REFPT | awk '{print $1-1}'`

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap F7/a/$BASE\_treeview.png \
	-l "-1" -m "0" -r "+1" -w 2 -f 18 \
	-x "Distance from TSS (kb)" -y "${NSITES} CoPRO determined TSSs sorted by CpG island length" \
	-o F7/a/$BASE\_treeview.svg

# Copy same file into supplement
[ -d S1 ] || mkdir S1
cp F7/a/$BASE\_treeview.svg S1/$BASE\_treeview.svg

# ===============================================================================================================================

BAMFILE=$WRK/../data/BAM/MNase-seq_ENCODE_merge_hg38.bam
BAM=`basename $BAMFILE ".bam"`

echo "Pileup MNase-seq across TSS RefPT w/ 80bp shift"
BASE=$BAM\_$BED

# Pileup SE MNase data (shift 80bp)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $REFPT $BAMFILE \
	-5 -1 --shift 80 --combined --cpu $THREADS -z \
	-o F7/a/$BASE\_combined.out -M F7/a/$BASE

# Two-color heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --black -p 0.95 F7/a/$BASE\_combined.cdt.gz -o F7/a/$BASE\_treeview.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap F7/a/$BASE\_treeview.png \
	-l "-1" -m "0" -r "+1" -w 2 -f 18 \
	-x "Distance from TSS (kb)" -y "${NSITES} CoPRO determined TSSs sorted by CpG island length" \
	-o F7/a/$BASE\_treeview.svg

# ===============================================================================================================================

echo "Copy heatmaps over from X_Bulk_Processing/Library"

MNASE_BAMFILE=../data/BAM/MNase-seq_ENCODE_merge_hg38.bam
BNASE_BAMFILE=../data/BAM/BNase-seq_50U-10min_merge_hg38.bam

DATAFILE=F7/a/violin_data.txt
[ -f $DATAFILE ] && rm $DATAFILE

for GROUP in "NoOverlap" "Overlap";
do
	BEDFILE=../data/RefPT-Krebs/PlusOneDyad_SORT-DistToExpressedTSS_GROUP-CpGIsland-$GROUP.bed

	# Tag pileup BNase-seq and MNase-seq
	java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 --combined --shift 80 $BEDFILE $MNASE_BAMFILE -M F7/a/MNase_$GROUP\_5read1-SHIFT80
	java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -p --combined            $BEDFILE $BNASE_BAMFILE -M F7/a/BNase_$GROUP\_midpoint
	# java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 --combined --shift 80 $BEDFILE $BNASE_BAMFILE -M F7/a/BNase_$GROUP\_5read1-SHIFT80

	# Aggregate data into violin data file
	java -jar $SCRIPTMANAGER read-analysis aggregate-data -m F7/a/MNase_$GROUP\_5read1-SHIFT80_combined.cdt -o F7/a/BNase_$GROUP.out
	java -jar $SCRIPTMANAGER read-analysis aggregate-data -m F7/a/BNase_$GROUP\_midpoint_combined.cdt       -o F7/a/MNase_$GROUP.out

	# Append aggregated data into a merged violin data file
	sed '1d' F7/a/BNase_$GROUP.out | awk -v GROUP=$GROUP 'BEGIN{OFS="\t";FS="\t"}{print $2,"BNase-"GROUP}' >> $DATAFILE
	sed '1d' F7/a/MNase_$GROUP.out | awk -v GROUP=$GROUP 'BEGIN{OFS="\t";FS="\t"}{print $2,"MNase-"GROUP}' >> $DATAFILE
done

# Violins
python $VIOLIN -i $DATAFILE -o F7/a/violin_data.png \
	--xlabel "+1 dyad group" --ylabel "Tag Occupancy" \
	--width 8 --height 4\
	--title "BNase v MNase in CpGIsland overlap v non-overlap regions"