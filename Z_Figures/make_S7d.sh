#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=14gb
#SBATCH -t 00:15:00
#SBATCH -A open
#SBATCH -o test.out
#SBATCH -e test.err

# Copy over composite data for S7/d from Library

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
[ -d S7/d ] || mkdir S7/d

# ===============================================================================================================================

# Composites
BED=PlusOneDyad_SORT-Expression_2000bp
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2AZ_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H2AZ.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H2B_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H2B.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H3_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H3.out
cp $LIBRARY/$BED/Composites/BNase-ChIP_H4_merge_hg38_$BED\_5read1-MIN128-MAX164.out S7/d/H4.out
cp $LIBRARY/$BED/Composites/BNase-seq_50U-10min_merge_hg38_$BED\_midpoint-MIN128-MAX164_combined.out S7/d/BNase.out
