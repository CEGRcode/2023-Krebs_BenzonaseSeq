#!/bin/bash

# Before digging into alignment and pre-processing scripts, simply copy existing
# merged BAM/BED files from Jordan into standard file naming system to work on
# building downstream figures.

#====Copy BED====
DIR=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/02_NFR_output_231128/
cp $DIR/K562_trueNFR.bed ../data/RefPT-Krebs/NFR_K562.bed
DIR=/storage/group/bfp2/default/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig1_atTSS_CpGsort/bedfiles/
cp $DIR/UCSCgb_hg19_CpGislands_230426.bed ../data/RefPT-Other/CpG_Islands.bed

#====Copy BAM====

BAMDIR=../data/BAM

# CoPRO
MERGED=/storage/group/bfp2/default/wkl2-WillLai/MEP_Project/01_BAM/CoPRO/
cp $MERGED/CoPRO_K562_MERGE.bam $BAMDIR/CoPRO_-_merge_hg19.bam
cp $MERGED/CoPRO_K562_Capped_SORT.bam $BAMDIR/CoPRO_Capped_merge_hg19.bam
cp $MERGED/CoPRO_K562_Uncapped_SORT.bam $BAMDIR/CoPRO_Uncapped_merge_hg19.bam

# XO
MERGED=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/
cp $MERGED/19354_19355_NoBenz_10sonicCycles_XO_PolII_master.bam $BAMDIR/ChIP-exo_Pol2_merge_hg19.bam

# BX
MERGED=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/
cp $MERGED/28452_28460_Benz_0sonicCycles_BX_H2A_master.bam $BAMDIR/BNase-ChIP_H2A_merge_hg19.bam
cp $MERGED/28453_28794_28461_28799_Benz_0sonicCycles_BX_H2A_Z_master.bam $BAMDIR/BNase-ChIP_H2AZ_merge_hg19.bam
cp $MERGED/28454_28795_28462_28800_Benz_0sonicCycles_BX_H2B_master.bam $BAMDIR/BNase-ChIP_H2B_merge_hg19.bam
cp $MERGED/25860_25868_25964_25971_28804_28808_Benz_0sonicCycles_BX_H3_master.bam $BAMDIR/BNase-ChIP_H3_merge_hg19.bam
cp $MERGED/28455_28463_Benz_0sonicCycles_BX_H3K4me1_master.bam $BAMDIR/BNase-ChIP_H3K4me1_merge_hg19.bam
cp $MERGED/25861_25869_25965_25972_28805_28809_Benz_0sonicCycles_BX_H3K4me3_master.bam $BAMDIR/BNase-ChIP_H3K4me3_merge_hg19.bam
cp $MERGED/25862_25870_25966_25973_Benz_0sonicCycles_BX_H3K9Ac_master.bam $BAMDIR/BNase-ChIP_H3K9ac_merge_hg19.bam
cp $MERGED/25858_25866_25962_25969_28806_28810_Benz_0sonicCycles_BX_H3K27Ac_master.bam $BAMDIR/BNase-ChIP_H3K27ac_merge_hg19.bam
cp $MERGED/28456_28796_28464_28801_Benz_0sonicCycles_BX_H3K27me3_master.bam $BAMDIR/BNase-ChIP_H3K27me3_merge_hg19.bam
cp $MERGED/25863_25871_25967_25974_Benz_0sonicCycles_BX_H3K36me3_master.bam $BAMDIR/BNase-ChIP_H3K36me3_merge_hg19.bam
cp $MERGED/28457_28797_28465_28802_Benz_0sonicCycles_BX_H4_master.bam $BAMDIR/BNase-ChIP_H4_merge_hg19.bam

# BI
BZNASE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/
cp $BZNASE/K562_benzonase-seq_master.bam $BAMDIR/BNase-seq_50U_merge_hg19.bam

# DNase sample & MNase titration
MNASE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_MNase_DNase/final_files/
cp $MNASE/SRR16815400_master.bam $BAMDIR/DNase-seq_-_rep1_hg19.bam
cp $MNASE/SRR3211679_master.bam $BAMDIR/MNase-seq_5U_rep1_hg19.bam
cp $MNASE/SRR3211680_master.bam $BAMDIR/MNase-seq_21U_rep1_hg19.bam
cp $MNASE/SRR3211681_master.bam $BAMDIR/MNase-seq_79U_rep1_hg19.bam
cp $MNASE/SRR3211682_master.bam $BAMDIR/MNase-seq_304U_rep1_hg19.bam

# Merged single-end deep MNase-seq
#??????

# MNase-ChIP (H3K4me3)
MCHIP=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240215_MNase_data/final_merged_files/
cp $MCHIP/SRR6010180_SRR6010175_SRR6010177_SRR7441419_SRR7441420_dedup_MERGE.bam $BAMDIR/MNase-ChIP_H3K4me3_merge_hg19.bam

# CUTandRUN
CUTRUN=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN
cp $CUTRUN/H2AZ/H2AZ_4DN_CUTRUN.bam $BAMDIR/CUTandRUN_H2AZ_merge_hg19.bam
cp $CUTRUN/H3K27ac/H3K27ac_4DN_CUTRUN.bam $BAMDIR/CUTandRUN_H3K27ac_merge_hg19.bam
cp $CUTRUN/H3K27me3/H3K27me3_4DN_CUTRUN.bam $BAMDIR/CUTandRUN_H3K27me3_merge_hg19.bam
cp $CUTRUN/H3K4me1/H3K4me1_4DN_CUTRUN.bam $BAMDIR/CUTandRUN_H3K4me1_merge_hg19.bam
cp $CUTRUN/H3K4me3/H3K4me3_4DN_CUTRUN.bam $BAMDIR/CUTandRUN_H3K4me3_merge_hg19.bam
cp $CUTRUN/IgG/IgG_4DN_CUTRUN.bam $BAMDIR/CUTandRUN_IgG_merge_hg19.bam

for BAMFILE in $BAMDIR/*.bam;
do
	samtools index $BAMFILE
done
