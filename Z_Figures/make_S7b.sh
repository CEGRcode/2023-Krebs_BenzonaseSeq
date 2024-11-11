#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o test.out
#SBATCH -e test.err

# Copy over composite data for S7/b from Library

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
[ -d S7/b ] || mkdir S7/b

# ===============================================================================================================================

# Heatmaps
BED=PlusOneDyad_SORT-Expression_2000bp
cp $LIBRARY/$BED/SVG/BNase-ChIP_H2A_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H2B_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K36me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b
cp $LIBRARY/$BED/SVG/BNase-ChIP_H4_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/b

# Composites
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2A_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H2A.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H2AZ.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2B_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H2B.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H3K4me1.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H3K27ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H3K27me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K36me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H3K36me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H4_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/b/H4.out
