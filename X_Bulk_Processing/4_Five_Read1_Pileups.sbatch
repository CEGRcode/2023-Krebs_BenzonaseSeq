#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/4_Five_Read1_Pileups.log.out-%a
#SBATCH -e logs/4_Five_Read1_Pileups.log.err-%a
#SBATCH --array 1-18

# Pileup 5' R1 (heatmaps) and composites for a custom combination of BAM x BED
# files with RGB info and custom insert size filters.
# Primarily ChIP data with some CUT&RUN.

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
METADATA=Five_Read1_Pileups.txt
THREADS=4
###

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
OUTDIR=$WRK/Library

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Determine RGB hexstring for the current job array index
RGB=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`

# Determine RGB hexstring for the current job array index
MIN=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
MAX=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Determine BAM file for the current job array index
BAMFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $4}'`
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index
BEDFILE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $5}'`
BED=`basename $BEDFILE ".bed"`

# Count sites
NSITES=`wc -l $BEDFILE | awk '{print $1}'`

# Set window labels
RLABEL="R"
LLABEL="L"
WINDOW=`echo $BED | rev | cut -d"_" -f1 | rev`
[ $WINDOW == "2000bp" ] && LLABEL="-1" && RLABEL="+1"
[ $WINDOW == "1000bp" ] && LLABEL="-500" && RLABEL="+500"
[ $WINDOW == "200bp" ] && LLABEL="-100" && RLABEL="+100"

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} "x" ${BAMFILE} "

DIR=$OUTDIR/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

# ===============================================================================================================================

echo "Run pileup (5' end, Read 1, combined)"
BASE=$BAM\_$BED\_5read1

[[ $MIN =~ "NaN" ]] || BASE=$BASE\-MIN$MIN
[[ $MIN =~ "NaN" ]] || PARAMS="-n $MIN"
[[ $MAX =~ "NaN" ]] || BASE=$BASE\-MAX$MAX
[[ $MAX =~ "NaN" ]] || PARAMS=$PARAMS" -x $MAX"

# Pileup (5read1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu $THREADS -5 -1 -p $PARAMS --combined -z -o $DIR/Composites/$BASE\_combined.out -M $DIR/CDT/$BASE
# No scaling

# Two-color heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -p .95 --color $RGB $DIR/CDT/$BASE\_combined.cdt.gz -o $DIR/PNG/$BASE\_combined_treeview.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_combined_treeview.png \
	-l $LLABEL -m "0" -r $RLABEL -w 2 -f 18 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/$BASE\_combined_treeview_label.svg

# ===============================================================================================================================

echo "Run pileup (5' end, Read 1, seperate)"
BASE=$BAM\_$BED\_5read1

[[ $MIN =~ "NaN" ]] || BASE=$BASE\-MIN$MIN || PARAMS="-n $MIN"
[[ $MAX =~ "NaN" ]] || BASE=$BASE\-MAX$MAX || PARAMS=$PARAMS" -x $MAX"

# Pileup (5read1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu $THREADS -5 -1 -p $PARAMS -z -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE
# No scaling

# Two-color heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --blue -p .95 $DIR/CDT/$BASE\_sense.cdt.gz -o $DIR/PNG/Strand/$BASE\_sense_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap --red  -p .95 $DIR/CDT/$BASE\_anti.cdt.gz  -o $DIR/PNG/Strand/$BASE\_anti_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/$BASE\_sense_treeview.png $DIR/PNG/Strand/$BASE\_anti_treeview.png -o $DIR/PNG/$BASE\_merge_treeview.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_merge_treeview.png \
	-l $LLABEL -m "0" -r $RLABEL -w 2 -f 18 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/$BASE\_merge_treeview_label.svg
