#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o test.out
#SBATCH -e test.err

# Copy over heatmap and composite data for S7/f from Library

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/Z_Figures
THREADS=4
###

set -exo

# Inputs and outputs
LIBRARY=$WRK/../X_Bulk_Processing/Library

[ -d S7 ] || mkdir S7
[ -d S7/f ] || mkdir S7/f

# ===============================================================================================================================

# Heatmaps
BED=PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp
cp $LIBRARY/$BED/SVG/CoPRO_Capped_merge_hg38_$BED\_5read2_merge_treeview_label.svg S7/f
cp $LIBRARY/$BED/SVG/CoPRO_Capped_merge_hg38_$BED\_5read1_merge_treeview_label.svg S7/f
cp $LIBRARY/$BED/SVG/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_combined_treeview_label.svg S7/f
cp $LIBRARY/$BED/SVG/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint_combined_treeview_label.svg S7/f

# Composites
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read2.out S7/f
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read1.out S7/f
cp $LIBRARY/$BED/Composites/ChIP-exo_Pol2_merge_hg38_$BED\_5read1_combined.out S7/f
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint_combined.out S7/f