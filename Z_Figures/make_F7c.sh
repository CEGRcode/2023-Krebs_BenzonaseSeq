#!/bin/bash

# Copy over composite data for F7/c from Library

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/Z_Figures
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/Z_Figures
THREADS=4
###

set -exo

# Inputs and outputs
LIBRARY=$WRK/../X_Bulk_Processing/Library

[ -d F7 ] || mkdir F7
[ -d F7/c ] || mkdir F7/c

# ===============================================================================================================================

BED=PlusOneDyad_SORT-Expression_2000bp

# Composites
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2A_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/c/BNase-ChIP-H2A.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2AZ_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/c/BNase-ChIP-H2AZ.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/c/BNase-ChIP-H3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/c/BNase-ChIP-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/c/BNase-ChIP-H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/c/BNase-ChIP-H3K27ac.out
cp $LIBRARY/$BED/Composites/MNase-ChIP_H3K4me3_merge_hg38_$BED\_midpoint-MAX80_combined.out F7/c/MNase-ChIP-H3K4me3.out
