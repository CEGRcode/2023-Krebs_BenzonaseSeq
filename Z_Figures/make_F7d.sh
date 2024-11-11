#!/bin/bash

# Copy over composite data for F7/d from Library

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/Z_Figures
THREADS=4
###

# Dependencies

set -exo

# Inputs and outputs
LIBRARY=$WRK/../X_Bulk_Processing/Library

# Script shortcuts
SCRIPTMANAGER=$WRK/../bin/ScriptManager-v0.15.jar

[ -d F7 ] || mkdir F7
[ -d F7/d ] || mkdir F7/d

# ===============================================================================================================================

# Heatmaps
BED=PlusOneDyad_SORT-pHN-dHN_400bp
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_combined_treeview_label.svg F7/d
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MAX80_combined_treeview_label.svg F7/d
cp $LIBRARY/$BED/SVG/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MAX80_combined_treeview_label.svg F7/d
java -jar -Djava.awt.headless=true $SCRIPTMANAGER figure-generation label-heatmap $LIBRARY/$BED/PNG/Strand/CoPRO_Capped_merge_hg38_$BED\_5read1_anti_treeview.png \
	-l -200 -m 0 -r +200 -w 2 -f 18 \
	-x $BED -y "$BED occurences (NSITES sites)" \
	-o F7/d/CoPRO_Capped_merge_hg38_$BED\_5read1_anti_treeview_label.svg

# Composites
BED=PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500_400bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/d/Top-H3K4me3.out
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read1.out F7/d/Top-CoPRO.out

BED=PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500_400bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/d/Bottom-H3K4me3.out
cp $LIBRARY/$BED/Composites/CoPRO_Capped_merge_hg38_$BED\_5read1.out F7/d/Bottom-CoPRO.out