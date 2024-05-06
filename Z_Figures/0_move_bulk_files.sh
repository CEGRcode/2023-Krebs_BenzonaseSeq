#!/bin/bash

# Organize data from X_Bulk_Processing into Z_Figures

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_Benzonase-seq/Z_Figures
###

LIBRARY=$WRK/../X_Bulk_Processing/Library

# ===============================================================================================================================

[ -d F1 ] || mkdir F1

[ -d F1/b ] || mkdir F1/b
# CpG heatmap - see custom script
# MNase heatmap - see custom script
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/BNase-seq_50U_merge_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg F1/b/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/DNase-seq_-_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg F1/b/
