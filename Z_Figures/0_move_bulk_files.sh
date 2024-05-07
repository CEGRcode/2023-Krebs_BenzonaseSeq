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

[ -d F4 ] || mkdir F4

[ -d F4/a ] || mkdir F4/a
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H2A_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H2AZ_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H2B_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H3K4me1_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H3K9ac_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H3K27ac_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H3K27me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H3K36me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/BNase-ChIP_H4_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/a/

[ -d F4/b ] || mkdir F4/b
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H2A_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H2AZ_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H2B_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3K4me1_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3K9ac_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3K27ac_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3K27me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3K36me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H4_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/b/

# Rename for plotter
for FILE in F4/b/*.out;
do
	BASE=`basename $FILE`
	TARGET=`echo $BASE | cut -d"_" -f2`
	mv $FILE F4/b/$TARGET\_$BASE
done

[ -d F4/c ] || mkdir F4/c
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-seq_50U_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint_combined.out F4/c/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H2AZ_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_5read1-MIN128-MAX164.out F4/c/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H2B_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_5read1-MIN128-MAX164.out F4/c/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_5read1-MIN128-MAX164.out F4/c/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/BNase-ChIP_H4_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_5read1-MIN128-MAX164.out F4/c/

# Rename for plotter
for FILE in F4/c/*.out;
do
	BASE=`basename $FILE`
	TARGET=`echo $BASE | cut -d"_" -f2`
	mv $FILE F4/c/$TARGET\_$BASE
done

[ -d F4/d ] || mkdir F4/d
# Violin plots of occupancy - see custom script

[ -d F4/e ] || mkdir F4/e
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/CUTandRUN_H2AZ_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/e/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/CUTandRUN_H3K4me1_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/e/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/CUTandRUN_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/e/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/CUTandRUN_H3K27ac_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/e/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/CUTandRUN_H3K27me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/e/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/SVG/CUTandRUN_IgG_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined_treeview_label.svg F4/e/

[ -d F4/f ] || mkdir F4/f
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CUTandRUN_H2AZ_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/f/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CUTandRUN_H3K4me1_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/f/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CUTandRUN_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/f/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CUTandRUN_H3K27ac_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/f/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CUTandRUN_H3K27me3_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/f/
cp $LIBRARY/PlusOneDyad_SORT-Expression_2000bp/Composites/CUTandRUN_IgG_merge_hg19_PlusOneDyad_SORT-Expression_2000bp_midpoint-MIN128-MAX164_combined.out F4/f/


# ===============================================================================================================================

[ -d F5 ] || mkdir F5

[ -d F5/a ] || mkdir F5/a
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/SVG/BNase-ChIP_H2A_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined_treeview_label.svg F5/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/SVG/BNase-ChIP_H2AZ_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined_treeview_label.svg F5/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/SVG/BNase-ChIP_H3_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined_treeview_label.svg F5/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/SVG/BNase-ChIP_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined_treeview_label.svg F5/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/SVG/BNase-ChIP_H3K9ac_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined_treeview_label.svg F5/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/SVG/BNase-ChIP_H3K27ac_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined_treeview_label.svg F5/a/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/SVG/MNase-ChIP_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined_treeview_label.svg F5/a/

[ -d F5/b ] || mkdir F5/b
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/Composites/BNase-ChIP_H2A_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined.out F5/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/Composites/BNase-ChIP_H2AZ_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined.out F5/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/Composites/BNase-ChIP_H3_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined.out F5/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/Composites/BNase-ChIP_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined.out F5/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/Composites/BNase-ChIP_H3K9ac_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined.out F5/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/Composites/BNase-ChIP_H3K27ac_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined.out F5/b/
cp $LIBRARY/PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp/Composites/MNase-ChIP_H3K4me3_merge_hg19_PlusOneDyad_SORT-Expression_GROUP-Nuc-Dyad_2000bp_midpoint-MIN92-MAX127_combined.out F5/b/


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
