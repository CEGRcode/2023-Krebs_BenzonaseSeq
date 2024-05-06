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

[ -d F3 ] || mkdir F3

[ -d F3/a ] || mkdir F3/a
cp $LIBRARY/PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp/SVG/CoPRO_-_merge_hg19_PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp_5read2_merge_treeview_label.svg F3/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp/SVG/CoPRO_-_merge_hg19_PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp_5read1_merge_treeview_label.svg F3/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp/SVG/ChIP-exo_Pol2_merge_hg19_PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp_5read1_combined_treeview_label.svg F3/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp/SVG/BNase-seq_50U_merge_hg19_PlusOneDyad_SORT-Expression_WithUnexpressed_2000bp_midpoint_combined_treeview_label.svg F3/a/

[ -d F3/b ] || mkdir F3/b
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CoPRO_-_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_5read2.out F3/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CoPRO_-_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_5read1.out F3/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/ChIP-exo_Pol2_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_5read1_combined.out F3/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-seq_50U_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint_combined.out F3/b/

# ===============================================================================================================================

[ -d S1 ] || mkdir S1
# CpG heatmap - see custom script
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_5U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_21U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_79U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/MNase-seq_304U_rep1_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/
cp $LIBRARY/TSS_GROUP-Expressed_SORT-CpG_2000bp/SVG/BNase-seq_50U_merge_hg19_TSS_GROUP-Expressed_SORT-CpG_2000bp_midpoint_combined_treeview_label.svg S1/

# ===============================================================================================================================

[ -d S2 ] || mkdir S2
# BNase-seq titration insert size histograms - see custom script
