#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o test.out
#SBATCH -e test.err

# Copy over composite data for S7/g from Library

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
[ -d S7/g ] || mkdir S7/g

# ===============================================================================================================================

# Composites
BED=PlusOneDyad_SORT-pHN-dHN_GROUP-TOP-2500_1000bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_5read1-MAX80.out S7/g/TOP-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_5read1-MAX80.out S7/g/TOP-H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_5read1-MAX80.out S7/g/TOP-H3K27ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_5read1-MAX80.out S7/g/TOP-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint-MAX80_combined.out S7/g/TOP-BNase.out

BED=PlusOneDyad_SORT-pHN-dHN_GROUP-BOTTOM-2500_1000bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K9ac_merge_hg38_$BED\_5read1-MAX80.out S7/g/BOTTOM-H3K9ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K27ac_merge_hg38_$BED\_5read1-MAX80.out S7/g/BOTTOM-H3K27ac.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3K4me3_merge_hg38_$BED\_5read1-MAX80.out S7/g/BOTTOM-H3K4me3.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint-MAX80_combined.out S7/g/BOTTOM-BNase.out
