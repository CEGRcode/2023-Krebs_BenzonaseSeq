#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o test.out
#SBATCH -e test.err

# Copy over heatmap and composite data for S7/c from Library

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
[ -d S7/c ] || mkdir S7/c

# ===============================================================================================================================

# Heatmaps
BED=PlusOneDyad_SORT-Expression_2000bp
cp $LIBRARY/$BED/SVG/CUTRUN_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/c
cp $LIBRARY/$BED/SVG/CUTRUN_IgG_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined_treeview_label.svg S7/c

# Composites
cp $LIBRARY/$BED/Composites/CUTRUN_H2AZ_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K4me1_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K4me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K27ac_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_H3K27me3_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/c
cp $LIBRARY/$BED/Composites/CUTRUN_IgG_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/c
