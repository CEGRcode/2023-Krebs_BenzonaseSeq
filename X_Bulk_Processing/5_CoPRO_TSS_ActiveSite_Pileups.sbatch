#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=24gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o logs/5_CoPRO_TSS_ActiveSite_Pileups.log.out-%a
#SBATCH -e logs/5_CoPRO_TSS_ActiveSite_Pileups.log.err-%a
#SBATCH --array 1-5

# Pileup 5' R1 and R2 (heatmaps) and composites for a list of BED files using
# the CoPRO BAM file and TSS/Active site heatmap colors.

# Scaling? Color threshold?

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
THREADS=4
###

LIST=(
	"../data/RefPT-Krebs/400bp/PlusOneDyad_SORT-pHN-dHN_400bp.bed"
	"../data/RefPT-Krebs/400bp/PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500_400bp.bed"
	"../data/RefPT-Krebs/400bp/PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500_400bp.bed"
	"../data/RefPT-Krebs/2000bp/PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp.bed"
	"../data/RefPT-Krebs/2000bp/PlusOneDyad_SORT-Expression_2000bp.bed"
)

# Dependencies
# - java
# - samtools

set -exo
module load samtools

# Fill in placeholder constants with your directories
BAMDIR=$WRK/../data/BAM
OUTDIR=$WRK/Library
BAMFILE=$BAMDIR/CoPRO_Capped_merge_hg38.bam

# Setup ScriptManager for job array
ORIGINAL_SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15-$SLURM_ARRAY_TASK_ID.jar
cp $ORIGINAL_SCRIPTMANAGER $SCRIPTMANAGER

# Script shortcuts
COMPOSITE=$WRK/../bin/sum_Col_CDT.pl

# Set up output directories
[ -d logs ] || mkdir logs
[ -d $OUTDIR ] || mkdir $OUTDIR

# Parse BAM filename and index if needed (CoPRO)
BAM=`basename $BAMFILE ".bam"`
[ -f $BAMFILE.bai ] || samtools index $BAMFILE

# Determine BED file for the current job array index (adjust for 0-indexed bash array)
INDEX=$(($SLURM_ARRAY_TASK_ID-1))
BEDFILE=${LIST[$INDEX]}
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
[ $WINDOW == "400bp" ] && LLABEL="-200" && RLABEL="+200"

echo "(${SLURM_ARRAY_TASK_ID}) ${BEDFILE} "x" ${BAMFILE} "

DIR=$OUTDIR/$BED
[ -d $DIR ] || mkdir $DIR
[[ -d $DIR/CDT ]] || mkdir $DIR/CDT
[[ -d $DIR/Composites ]] || mkdir $DIR/Composites
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

# ===============================================================================================================================

echo "Run CoPRO active site pileup (5' end, Read 1)"
BASE=$BAM\_$BED\_5read1

# Pileup (5read1)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu $THREADS -5 -1 -p -z -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE
# No scaling

# Two-color heatmap -a 0.4
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -a 1 -c 7F00FF $RGB $DIR/CDT/$BASE\_sense.cdt.gz -o $DIR/PNG/Strand/$BASE\_sense_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -a 1 -c 0099FF $RGB $DIR/CDT/$BASE\_anti.cdt.gz -o $DIR/PNG/Strand/$BASE\_anti_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/$BASE\_sense_treeview.png $DIR/PNG/Strand/$BASE\_anti_treeview.png -o $DIR/PNG/$BASE\_merge_treeview.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_merge_treeview.png \
	-l $LLABEL -m "0" -r $RLABEL -w 2 -f 18 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/$BASE\_merge_treeview_label.svg

# ===============================================================================================================================

echo "Run CoPRO TSS pileup (5' end, Read 2)"
BASE=$BAM\_$BED\_5read2

# Pileup (5read2)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu $THREADS -5 -2 -p -z -o $DIR/Composites/$BASE.out -M $DIR/CDT/$BASE
# No scaling

# Two-color heatmap (?) -a 0.6
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -a 1 -c FF0000 $RGB $DIR/CDT/$BASE\_sense.cdt.gz -o $DIR/PNG/Strand/$BASE\_sense_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation heatmap -a 1 -c 0000FF $RGB $DIR/CDT/$BASE\_anti.cdt.gz -o $DIR/PNG/Strand/$BASE\_anti_treeview.png
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation merge-heatmap $DIR/PNG/Strand/$BASE\_sense_treeview.png $DIR/PNG/Strand/$BASE\_anti_treeview.png -o $DIR/PNG/$BASE\_merge_treeview.png

# Label heatmap
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $DIR/PNG/$BASE\_merge_treeview.png \
	-l $LLABEL -m "0" -r $RLABEL -w 2 -f 18 \
	-x $BED -y "$BED occurences (${NSITES} sites)" \
	-o $DIR/SVG/$BASE\_merge_treeview_label.svg