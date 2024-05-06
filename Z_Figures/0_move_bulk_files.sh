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

[ -d F1/c ] || mkdir F1/c
# Insert (BI) - see custom script

[ -d F1/e ] || mkdir F1/e
# Insert (MNase 21U) - see custom script
# Insert (MNase 304U) - see custom script

[ -d F1/f ] || mkdir F1/f
# Insert (DNase) - see custom script

# ===============================================================================================================================

[ -d S1 ] || mkdir S1
# CpG heatmap - see custom script
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_5U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_21U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_79U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_304U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/BNase-seq_50U_merge_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
