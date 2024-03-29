#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:40:00
#SBATCH -A open
#SBATCH -o logs/2_MotifAnalyses.log.out-%a
#SBATCH -e logs/2_MotifAnalyses.log.err-%a
#SBATCH --array 1-13

# Loop through each BEDFILE in data/RefPT-Motif/1000bp and perform motif analyses:
# - Four color plots (QC)
# - BI pileups
# - DNA shape
# - Dinucleotide

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/X_Bulk_Processing
METADATA=BI_Pileups.txt
THREADS=4
###

# Dependencies
# - java
# - opencv
# - python

set -exo
module load anaconda
conda activate bx

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.14.jar
RESIZE=../bin/resize_png.py
DINUCLEOTIDE=../bin/dinucleotide_CDT_from_FASTA.py
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
[[ -d $DIR/FASTA ]] || mkdir $DIR/FASTA
[[ -d $DIR/FourColor ]] || mkdir $DIR/FourColor
[[ -d $DIR/PNG/Strand ]] || mkdir -p $DIR/PNG/Strand
[[ -d $DIR/SVG ]] || mkdir $DIR/SVG

# ===============================================================================================================================

echo "Make Four-color plot"
# Expand 32bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 32 $BEDFILE -o $DIR/FourColor/$BED\_32bp.bed
# Filter for random 50 sites and take off 1bp (since all motifs are odd-lengthed, we can take 1bp off to center 0bp label on motif center)
shuf $DIR/FourColor/$BED\_32bp.bed | head -n 50 \
	| awk '{OFS="\t"}{if($6=="-") print $1,$2+1,$3,$4,$5,$6; else print $1,$2,$3-1,$4,$5,$6}' \
	> $DIR/FourColor/$BED\_31bp.bed
# Extract sequence
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $DIR/FourColor/$BED\_31bp.bed -o $DIR/FourColor/$BED\_31bp.fa
# Generate four-color plot from sequence
java -jar $SCRIPTMANAGER figure-generation four-color $DIR/FourColor/$BED\_31bp.fa -o $DIR/FourColor/$BED\_31bp.png
# Resize PNG to standard px dimensions
python $RESIZE -r 600 -c 200 -i $DIR/FourColor/$BED\_31bp.png -o $DIR/FourColor/$BED\_31bp_RESIZE.png
# Add SVG label
java -jar $SCRIPTMANAGER figure-generation label-heatmap $DIR/FourColor/$BED\_31bp_RESIZE.png \
	-l "-15" -m "0" -r "+15" -w 1 -f 20 \
	-x $BED -y "$BED Motif occurences (${NSITES} sites)" \
	-o $DIR/FourColor/$BED\_31bp.svg

# ===============================================================================================================================

echo "Run Motif endo cuts BI pileup"
BASE=$BAM\_$BED\_5both

# Pileup (endo)
java -jar $SCRIPTMANAGER read-analysis tag-pileup $BEDFILE $BAMFILE --cpu $THREADS -5 -a -o $DIR/Composites/$BASE.out
# No scaling

# ===============================================================================================================================

echo "Run DNA Shape analysis"

# Run DNA Shape
java -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed --all $GENOME $BEDFILE -o $DIR/CDT/$BED\_DNAShape

# CDT to composite
perl $COMPOSITE $DIR/CDT/$BED\_DNAShape_HelT.cdt  $DIR/Composites/DNAShape-HelT_$BED.out
perl $COMPOSITE $DIR/CDT/$BED\_DNAShape_MGW.cdt  $DIR/Composites/DNAShape-MGW_$BED.out
perl $COMPOSITE $DIR/CDT/$BED\_DNAShape_PropT.cdt  $DIR/Composites/DNAShape-PropT_$BED.out
perl $COMPOSITE $DIR/CDT/$BED\_DNAShape_Roll.cdt  $DIR/Composites/DNAShape-Roll_$BED.out

# ===============================================================================================================================

echo "Run Dinucleotide analysis"

# Extract sequence
java -jar $SCRIPTMANAGER sequence-analysis fasta-extract $GENOME $BEDFILE -o $DIR/FASTA/$BED.fa

for SEQ in "AA-AT-TA-TT" "CC-CG-GC-GG" "AC-CA-TG-GT" "AG-GA-TC-CT";
do
	# Dinucleotide analysis
	python $DINUCLEOTIDE -i $DIR/FASTA/$BED.fa --seq $SEQ -o $DIR/CDT/Dinucleotide-$SEQ\_$BED.cdt

	# CDT to composite
	perl $COMPOSITE $DIR/CDT/Dinucleotide-$SEQ\_$BED.cdt  $DIR/Composites/Dinucleotide-$SEQ\_$BED.out
done
