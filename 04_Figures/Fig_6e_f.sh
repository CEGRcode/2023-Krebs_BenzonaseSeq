# purpose - make Figs. 4c showing H3K4me3 and CoPRO at ratio of proximal to distal halves Nucleosomes at +1 nucleosome. *this version is based on sort of midpoint of reads. **comp plots are based top and bottom 2500 sites. ***Code checked by JEK 240209.

# usage
# qq
#
# example
# purpose - qq

# usage
# qq
#
# example
#
# 'qq'

#set bedfiles with Plus1 (downstream) nucleosome positions based on non-redundant bedfile-based calls
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig4_5_new_subnucleosomes/Fig5_H3K4me3_HN_sort_PolII_rebuilt_CoPRO_v8_240304

#set bam library file to BI_rep1 **testing with subsampled master BAM file
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25861_25869_25965_25972_28805_28809_Benz_0sonicCycles_BX_H3K4me3_master.bam
BAM2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25862_25870_25966_25973_Benz_0sonicCycles_BX_H3K9Ac_master.bam
BAM3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25858_25866_25962_25969_28806_28810_Benz_0sonicCycles_BX_H3K27Ac_master.bam
BAM4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/CoPRO/CoPRO_K562_MERGE.bam
BAM5=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/K562_benzonase-seq_master.bam
BAM6=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240215_MNase_data/final_merged_files/SRR6010180_SRR6010175_SRR6010177_SRR7441419_SRR7441420_dedup_MERGE_master.bam

#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl
JOB_ROW=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig6_Subnucleosomes/job/sum_Row_CDT.pl

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24GB
#SBATCH --time=4:00:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
module load anaconda
conda activate bioinfo"

#set output file names
BEDFILE_shuffled=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
BEDFILE_shuffled_a=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled"}')
BEDFILE_2bp=$(echo $BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2bp.bed"}')
Plus1_proximal_SN=$(echo $BEDFILE_2bp| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_proximal_SN.bed"}')
Plus1_distal_SN=$(echo $BEDFILE_2bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_distal_SN.bed"}')
Plus1_proximal_SN_a=$(echo $Plus1_proximal_SN | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
Plus1_distal_SN_a=$(echo $Plus1_distal_SN | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal"}')
CDT1_gz=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal_combined.cdt.gz"}')
CDT1_unzipped=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_proximal_combined.cdt"}')
OUT2=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal_midpoint.out"}')
CDT2=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal"}')
CDT2_gz=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal_combined.cdt.gz"}')
CDT2_unzipped=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_distal_combined.cdt"}')
Target_proximal=$(echo "$BAM1a""_""$Plus1_proximal_SN_a" | awk -F. '{print $1"_RowCount.tab"}')
Target_distal=$(echo "$BAM1a""_""$Plus1_distal_SN_a" | awk -F. '{print $1"_RowCount.tab"}')
RATIO_BEDFILE=$(echo "$BEDFILE_shuffled_a" | awk -F. '{print $1"_ratio_proximal_to_distal.bed"}')
RATIO_BEDFILE_2000=$(echo $RATIO_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2000bp.bed"}')
RATIO_BEDFILE_3000=$(echo $RATIO_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_3000bp.bed"}')
RATIO_BEDFILE_3000_top2500=$(echo $RATIO_BEDFILE_3000| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_top2500.bed"}')
RATIO_BEDFILE_3000_bottom2500=$(echo $RATIO_BEDFILE_3000| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_bottom2500.bed"}')
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM4a=$(echo $BAM4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM5a=$(echo $BAM5 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM6a=$(echo $BAM6 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
SCALE3=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT"}')
SCALE3a=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
SCALE4=$(echo "$BAM2a" | awk -F. '{print $1"_read2_ForCDT"}')
SCALE4a=$(echo "$BAM2a" | awk -F. '{print $1"_read2_ForCDT_ScalingFactors.out"}')
SCALE5=$(echo "$BAM3a" | awk -F. '{print $1"_read1_ForCDT"}')
SCALE5a=$(echo "$BAM3a" | awk -F. '{print $1"_read1_ForCDT_ScalingFactors.out"}')
SCALE6=$(echo "$BAM4a" | awk -F. '{print $1"_ForCDT"}')
SCALE6a=$(echo "$BAM4a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
SCALE7=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT"}')
SCALE7a=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
SCALE8=$(echo "$BAM6a" | awk -F. '{print $1"_ForCDT"}')
SCALE8a=$(echo "$BAM6a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
OUT3=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT3=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1}')
CDT3b=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT3c=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT3_SCALED=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
PNG3=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG3=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT3a_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT3d_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite"}')
CDT3e_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT3f_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT3_SCALED_COMP_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT3_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT4=$(echo "$BAM2a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT4=$(echo "$BAM2a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1}')
CDT4b=$(echo "$BAM2a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT4c=$(echo "$BAM2a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT4_SCALED=$(echo "$BAM2a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
PNG4=$(echo "$BAM2a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG4=$(echo "$BAM2a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT4a_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT4d_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite"}')
CDT4e_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT4f_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT4_SCALED_COMP_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT4_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT5=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT5=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1}')
CDT5b=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT5c=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT5_SCALED=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
PNG5=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG5=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT5a_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT5d_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite"}')
CDT5e_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT5f_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT5_SCALED_COMP_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT5_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT6=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1.out"}')
CDT6=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1"}')
CDT6_sense_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT6_anti_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT6_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT6_SCALED_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1_anti_scaled.cdt"}')
PNG6_sense=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1_scaled_sense.png"}')
PNG6_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read1_scaled_anti.png"}')
SVG6=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled_read1_labeled.svg"}')
OUT6a_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read1_ForComposite.out"}')
CDT6d_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read1_ForComposite"}')
CDT6_sense_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read1_ForComposite_sense.cdt.gz"}')
CDT6_anti_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read1_ForComposite_anti.cdt.gz"}')
CDT6_anti_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read1_ForComposite_anti.cdt"}')
CDT6_SCALED_COMP_anti_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.cdt"}')
SCALED_OUT6_anti_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.tab"}')
OUT12a_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read2_ForComposite.out"}')
CDT12d_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read2_ForComposite"}')
CDT12_sense_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read2_ForComposite_sense.cdt.gz"}')
CDT12_anti_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read2_ForComposite_anti.cdt.gz"}')
CDT12_sense_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read2_ForComposite_sense.cdt"}')
CDT12_SCALED_COMP_sense_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read2_ForComposite_scaled_sense.cdt"}')
SCALED_OUT12_sense_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_read2_ForComposite_scaled_sense.tab"}')
OUT12=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2.out"}')
CDT12=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2"}')
CDT12_sense_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_sense.cdt.gz"}')
CDT12_anti_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_anti.cdt.gz"}')
CDT12_sense=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_sense.cdt"}')
CDT12_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_anti.cdt"}')
CDT12_SCALED_sense=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_sense_scaled.cdt"}')
CDT12_SCALED_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_anti_scaled.cdt"}')
PNG12_sense=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_scaled_sense.png"}')
PNG12_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_scaled_anti.png"}')
PNG12_merge=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_read2_scaled_merged.png"}')
SVG12=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled_read2_labeled.svg"}')
OUT7=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT7=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1}')
CDT7b=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT7c=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT7_SCALED=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
PNG7=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG7=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT8=$(echo "$BAM6a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT8=$(echo "$BAM6a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1}')
CDT8b=$(echo "$BAM6a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT8c=$(echo "$BAM6a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT8_SCALED=$(echo "$BAM6a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
PNG8=$(echo "$BAM6a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG8=$(echo "$BAM6a""_""$RATIO_BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT8a_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT8d_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite"}')
CDT8e_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT8f_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT8_SCALED_COMP_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT8_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_top2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT23a_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT23d_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite"}')
CDT23e_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT23f_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT23_SCALED_COMP_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT23_3000=$(echo "$BAM1a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT24a_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT24d_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite"}')
CDT24e_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT24f_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT24_SCALED_COMP_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT24_3000=$(echo "$BAM2a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT25a_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT25d_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite"}')
CDT25e_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT25f_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT25_SCALED_COMP_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT25_3000=$(echo "$BAM3a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT26a_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read1_ForComposite.out"}')
CDT26d_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read1_ForComposite"}')
CDT26_sense_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read1_ForComposite_sense.cdt.gz"}')
CDT26_anti_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read1_ForComposite_anti.cdt.gz"}')
CDT26_anti_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read1_ForComposite_anti.cdt"}')
CDT26_SCALED_COMP_anti_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.cdt"}')
SCALED_OUT26_anti_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.tab"}')
OUT32a_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read2_ForComposite.out"}')
CDT32d_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read2_ForComposite"}')
CDT32_sense_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read2_ForComposite_sense.cdt.gz"}')
CDT32_anti_3000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read2_ForComposite_anti.cdt.gz"}')
CDT32_sense_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read2_ForComposite_sense.cdt"}')
CDT32_SCALED_COMP_sense_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read2_ForComposite_scaled_sense.cdt"}')
SCALED_OUT32_sense_3000=$(echo "$BAM4a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_read2_ForComposite_scaled_sense.tab"}')
OUT28a_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT28d_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite"}')
CDT28e_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT28f_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT28_SCALED_COMP_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT28_3000=$(echo "$BAM6a""_""$RATIO_BEDFILE_3000_bottom2500" | awk -F. '{print $1"_ForComposite_scaled.tab"}')

sampleID=fig5_H3K4me3_prox_distal_ratio_CoPRO_rebuilt_v8_240304.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#randomize intitial bedfile so that there any no secondary sorts" >> $sampleID
echo "shuf $BEDFILE > $BEDFILE_shuffled" >> $sampleID
echo "#expand bedfiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2 $BEDFILE_shuffled -o=$BEDFILE_2bp" >> $sampleID
echo "#get proximal (upstream) or distal (downstream) half of nucleosome relative to +1 midpoint (dyad)" >> $sampleID
echo "#math is based on 80 bp SNs for shifts: proximal SN is ~2 bp upstream of dyad so 80 + 2 = 82 for start, -2 + -2 = -4 for end; distal SN is ~4 bp downstream of dyad so -2 + -4 = -6 for start, 80 + 4 = 84 for end" >> $sampleID
echo "bedtools slop -i $BEDFILE_2bp -g $HG19_GENOME -l 82 -r -4 -s > $Plus1_proximal_SN" >> $sampleID
echo "bedtools slop -i $BEDFILE_2bp -g $HG19_GENOME -l -6 -r 84 -s > $Plus1_distal_SN" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 --max-insert=80 $Plus1_proximal_SN $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 --max-insert=80 $Plus1_distal_SN $BAM1" >> $sampleID
echo "#unzip files" >> $sampleID
echo "gunzip -c $CDT1_gz > $CDT1_unzipped" >> $sampleID
echo "gunzip -c $CDT2_gz > $CDT2_unzipped" >> $sampleID
echo "#sum tag in each above CDT file" >> $sampleID
echo "perl $JOB_ROW $CDT1_unzipped $Target_proximal" >> $sampleID
echo "perl $JOB_ROW $CDT2_unzipped $Target_distal" >> $sampleID
echo "#only keep rows if all IDs match between initial bedfile, prox SN, and distal SN; then remove rows if 0 tags in numerator or denominator, then do division, sort by ratio and print bedfile" >> $sampleID
echo "paste $BEDFILE_shuffled $Target_proximal $Target_distal | awk '(\$4==\$7 && \$4==\$9 && \$8!=0 && \$10!=0){print \$0\"\t\"(\$8/\$10)}' | sort -k11,11rn | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $RATIO_BEDFILE" >> $sampleID
echo "paste $BEDFILE_shuffled $Target_proximal $Target_distal | awk '(\$4==\$7 && \$4==\$9 && \$8!=0 && \$10!=0){print \$0\"\t\"(\$8/\$10)}' | sort -k11,11rn > test_sort.bed" >> $sampleID
echo "#expand new bedfile by 2000 bp and 3000 bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2000 $RATIO_BEDFILE -o=$RATIO_BEDFILE_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=3000 $RATIO_BEDFILE -o=$RATIO_BEDFILE_3000" >> $sampleID
echo "#take top 2500 sites, middle 5000 sites, and bottom 2500 sites from above above bedfile. Bottom 2500 + ~180 rows that have PolII chunk -> take 2680 but leave last 180" >> $sampleID
echo "cat $RATIO_BEDFILE_3000 | head -2500 > $RATIO_BEDFILE_3000_top2500" >> $sampleID
echo "cat $RATIO_BEDFILE_3000 | tail -3073 | head -2500 > $RATIO_BEDFILE_3000_bottom2500" >> $sampleID
echo "#scale all libraries: options - total tag scaling -t" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE4 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE5 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE6 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE7 $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE8 $BAM6" >> $sampleID
echo "#do tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 --max-insert=80 $RATIO_BEDFILE_2000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 --max-insert=80 $RATIO_BEDFILE_2000 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 --max-insert=80 $RATIO_BEDFILE_2000 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT6 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6 $RATIO_BEDFILE_2000 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 -z --output-matrix=$CDT12 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT12 $RATIO_BEDFILE_2000 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT7 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT7 $RATIO_BEDFILE_2000 $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT8 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT8 --max-insert=80 $RATIO_BEDFILE_2000 $BAM6" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT3b > $OUTPUT/$CDT3c" >> $sampleID
echo "gunzip -c $CDT4b > $OUTPUT/$CDT4c" >> $sampleID
echo "gunzip -c $CDT5b > $OUTPUT/$CDT5c" >> $sampleID
echo "gunzip -c $CDT6_anti_gz > $CDT6_anti" >> $sampleID
echo "gunzip -c $CDT12_sense_gz > $CDT12_sense" >> $sampleID
echo "gunzip -c $CDT12_anti_gz > $CDT12_anti" >> $sampleID
echo "gunzip -c $CDT7b > $OUTPUT/$CDT7c" >> $sampleID
echo "gunzip -c $CDT8b > $OUTPUT/$CDT8c" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT3c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED --scaling-factor=\$(cat $SCALE4a | cut -f2 | tail -1 | awk '{print \$1}') $CDT4c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT5c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_anti --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT12_SCALED_sense --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT12_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT12_SCALED_anti --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT12_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT7_SCALED --scaling-factor=\$(cat $SCALE7a | cut -f2 | tail -1 | awk '{print \$1}') $CDT7c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT8_SCALED --scaling-factor=\$(cat $SCALE8a | cut -f2 | tail -1 | awk '{print \$1}') $CDT8c" >> $sampleID
echo "#make heatmap [default black colors used] **change threshold if necessary" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG3 -c ff9900 -p=0.95 $CDT3_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG4 -c 0068c9 -p=0.95 $CDT4_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG5 -c 9900ff -p=0.95 $CDT5_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG6_anti -c 0099ff -p=0.95 $CDT6_SCALED_anti" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG12_sense -c ff0000 -p=0.95 $CDT12_SCALED_sense" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG12_anti -c 0000ff -p=0.95 $CDT12_SCALED_anti" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG7 -c ff00ff -p=0.95 $CDT7_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG8 -c 000000 -p=0.95 $CDT8_SCALED" >> $sampleID
echo "#merge heatmaps" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation merge-heatmap -o=$PNG12_merge $PNG12_sense $PNG12_anti" >> $sampleID
echo "#label above heatmaps" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG3 --output=$SVG3 --width=2 --font-size=18 --left-label='-1.0' --mid-label=0 --right-label='+1.0' --x-label='Distance from +1 Dyad (kb)' --y-label='10,733 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG4 --output=$SVG4 --width=2 --font-size=18 --left-label='-1.0' --mid-label=0 --right-label='+1.0' --x-label='Distance from +1 Dyad (kb)' --y-label='10,733 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG5 --output=$SVG5 --width=2 --font-size=18 --left-label='-1.0' --mid-label=0 --right-label='+1.0' --x-label='Distance from +1 Dyad (kb)' --y-label='10,733 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG6_anti --output=$SVG6 --width=2 --font-size=18 --left-label='-1.0' --mid-label=0 --right-label='+1.0' --x-label='Distance from +1 Dyad (kb)' --y-label='10,733 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG12_merge --output=$SVG12 --width=2 --font-size=18 --left-label='-1.0' --mid-label=0 --right-label='+1.0' --x-label='Distance from +1 Dyad (kb)' --y-label='10,733 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG7 --output=$SVG7 --width=2 --font-size=18 --left-label='-1.0' --mid-label=0 --right-label='+1.0' --x-label='Distance from +1 Dyad (kb)' --y-label='10,733 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG8 --output=$SVG8 --width=2 --font-size=18 --left-label='-1.0' --mid-label=0 --right-label='+1.0' --x-label='Distance from +1 Dyad (kb)' --y-label='10,733 Nucleosome positions sorted by ratio'" >> $sampleID
echo "#prep files for plotter for top2500 sites" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3a_3000 --max-insert=80 $RATIO_BEDFILE_3000_top2500 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4a_3000 --max-insert=80 $RATIO_BEDFILE_3000_top2500 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5a_3000 --max-insert=80 $RATIO_BEDFILE_3000_top2500 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT6d_3000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6a_3000 $RATIO_BEDFILE_3000_top2500 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 -z --output-matrix=$CDT12d_3000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT12a_3000 $RATIO_BEDFILE_3000_top2500 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT8d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT8a_3000 --max-insert=80 $RATIO_BEDFILE_3000_top2500 $BAM6" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT3e_3000 > $OUTPUT/$CDT3f_3000" >> $sampleID
echo "gunzip -c $CDT4e_3000 > $OUTPUT/$CDT4f_3000" >> $sampleID
echo "gunzip -c $CDT5e_3000 > $OUTPUT/$CDT5f_3000" >> $sampleID
echo "gunzip -c $CDT6_anti_3000_gz > $CDT6_anti_3000" >> $sampleID
echo "gunzip -c $CDT12_sense_3000_gz > $CDT12_sense_3000" >> $sampleID
echo "gunzip -c $CDT8e_3000 > $OUTPUT/$CDT8f_3000" >> $sampleID
echo "#scale CDT files data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT3f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE4a | cut -f2 | tail -1 | awk '{print \$1}') $CDT4f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT5f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_COMP_anti_3000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_anti_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT12_SCALED_COMP_sense_3000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT12_sense_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT8_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE8a | cut -f2 | tail -1 | awk '{print \$1}') $CDT8f_3000" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_3000 $SCALED_OUT3_3000" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_3000 $SCALED_OUT4_3000" >> $sampleID
echo "perl $JOB $CDT5_SCALED_COMP_3000 $SCALED_OUT5_3000" >> $sampleID
echo "perl $JOB $CDT6_SCALED_COMP_anti_3000 $SCALED_OUT6_anti_3000" >> $sampleID
echo "perl $JOB $CDT12_SCALED_COMP_sense_3000 $SCALED_OUT12_sense_3000" >> $sampleID
echo "perl $JOB $CDT8_SCALED_COMP_3000 $SCALED_OUT8_3000" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT23d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT23a_3000 --max-insert=80 $RATIO_BEDFILE_3000_bottom2500 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT24d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT24a_3000 --max-insert=80 $RATIO_BEDFILE_3000_bottom2500 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT25d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT25a_3000 --max-insert=80 $RATIO_BEDFILE_3000_bottom2500 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT26d_3000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT26a_3000 $RATIO_BEDFILE_3000_bottom2500 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 -z --output-matrix=$CDT32d_3000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT32a_3000 $RATIO_BEDFILE_3000_bottom2500 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT28d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT28a_3000 --max-insert=80 $RATIO_BEDFILE_3000_bottom2500 $BAM6" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT23e_3000 > $OUTPUT/$CDT23f_3000" >> $sampleID
echo "gunzip -c $CDT24e_3000 > $OUTPUT/$CDT24f_3000" >> $sampleID
echo "gunzip -c $CDT25e_3000 > $OUTPUT/$CDT25f_3000" >> $sampleID
echo "gunzip -c $CDT26_anti_3000_gz > $CDT26_anti_3000" >> $sampleID
echo "gunzip -c $CDT32_sense_3000_gz > $CDT32_sense_3000" >> $sampleID
echo "gunzip -c $CDT28e_3000 > $OUTPUT/$CDT28f_3000" >> $sampleID
echo "#scale CDT files data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT23_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT23f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT24_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE4a | cut -f2 | tail -1 | awk '{print \$1}') $CDT24f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT25_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT25f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT26_SCALED_COMP_anti_3000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT26_anti_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT32_SCALED_COMP_sense_3000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT32_sense_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT28_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE8a | cut -f2 | tail -1 | awk '{print \$1}') $CDT28f_3000" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT23_SCALED_COMP_3000 $SCALED_OUT23_3000" >> $sampleID
echo "perl $JOB $CDT24_SCALED_COMP_3000 $SCALED_OUT24_3000" >> $sampleID
echo "perl $JOB $CDT25_SCALED_COMP_3000 $SCALED_OUT25_3000" >> $sampleID
echo "perl $JOB $CDT26_SCALED_COMP_anti_3000 $SCALED_OUT26_anti_3000" >> $sampleID
echo "perl $JOB $CDT32_SCALED_COMP_sense_3000 $SCALED_OUT32_sense_3000" >> $sampleID
echo "perl $JOB $CDT28_SCALED_COMP_3000 $SCALED_OUT28_3000" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $Plus1_proximal_SN" >> $sampleID
echo "rm $Plus1_distal_SN" >> $sampleID
echo "rm $Plus1_proximal_SN_a" >> $sampleID
echo "rm $Plus1_distal_SN_a" >> $sampleID
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1_gz" >> $sampleID
echo "rm $CDT1_unzipped" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2_gz" >> $sampleID
echo "rm $CDT2_unzipped" >> $sampleID
echo "rm $Target_proximal" >> $sampleID
echo "rm $Target_distal" >> $sampleID
echo "rm $RATIO_BEDFILE" >> $sampleID
echo "rm $RATIO_BEDFILE_2000" >> $sampleID
echo "rm $RATIO_BEDFILE_3000" >> $sampleID
echo "rm $RATIO_BEDFILE_3000_top2500" >> $sampleID
echo "rm $RATIO_BEDFILE_3000_bottom2500" >> $sampleID
echo "rm $SCALE3a" >> $sampleID
echo "rm $SCALE4a" >> $sampleID
echo "rm $SCALE5a" >> $sampleID
echo "rm $SCALE6a" >> $sampleID
echo "rm $SCALE7a" >> $sampleID
echo "rm $SCALE8a" >> $sampleID
echo "rm $OUT3" >> $sampleID
echo "rm $CDT3b" >> $sampleID
echo "rm $CDT3c" >> $sampleID
echo "rm $CDT3_SCALED" >> $sampleID
echo "rm $PNG3" >> $sampleID
echo "rm $OUT3a_3000" >> $sampleID
echo "rm $CDT3e_3000" >> $sampleID
echo "rm $CDT3f_3000" >> $sampleID
echo "rm $CDT3_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT4" >> $sampleID
echo "rm $CDT4b" >> $sampleID
echo "rm $CDT4c" >> $sampleID
echo "rm $CDT4_SCALED" >> $sampleID
echo "rm $PNG4" >> $sampleID
echo "rm $OUT4a_3000" >> $sampleID
echo "rm $CDT4e_3000" >> $sampleID
echo "rm $CDT4f_3000" >> $sampleID
echo "rm $CDT4_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT5" >> $sampleID
echo "rm $CDT5b" >> $sampleID
echo "rm $CDT5c" >> $sampleID
echo "rm $CDT5_SCALED" >> $sampleID
echo "rm $PNG5" >> $sampleID
echo "rm $OUT5a_3000" >> $sampleID
echo "rm $CDT5e_3000" >> $sampleID
echo "rm $CDT5f_3000" >> $sampleID
echo "rm $CDT5_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT6" >> $sampleID
echo "rm $CDT6_sense_gz" >> $sampleID
echo "rm $CDT6_anti_gz" >> $sampleID
echo "rm $CDT6_anti" >> $sampleID
echo "rm $CDT6_SCALED_anti" >> $sampleID
echo "rm $PNG6_sense" >> $sampleID
echo "rm $PNG6_anti" >> $sampleID
echo "rm $OUT6a_3000" >> $sampleID
echo "rm $CDT6_sense_3000_gz" >> $sampleID
echo "rm $CDT6_anti_3000_gz" >> $sampleID
echo "rm $CDT6_anti_3000" >> $sampleID
echo "rm $CDT6_SCALED_COMP_anti_3000" >> $sampleID
echo "rm $OUT12a_3000" >> $sampleID
echo "rm $CDT12_sense_3000_gz" >> $sampleID
echo "rm $CDT12_anti_3000_gz" >> $sampleID
echo "rm $CDT12_sense_3000" >> $sampleID
echo "rm $CDT12_SCALED_COMP_sense_3000" >> $sampleID
echo "rm $OUT12" >> $sampleID
echo "rm $CDT12_sense_gz" >> $sampleID
echo "rm $CDT12_anti_gz" >> $sampleID
echo "rm $CDT12_sense" >> $sampleID
echo "rm $CDT12_anti" >> $sampleID
echo "rm $CDT12_SCALED_sense" >> $sampleID
echo "rm $CDT12_SCALED_anti" >> $sampleID
echo "rm $PNG12_sense" >> $sampleID
echo "rm $PNG12_anti" >> $sampleID
echo "rm $PNG12_merge" >> $sampleID
echo "rm $OUT7" >> $sampleID
echo "rm $CDT7b" >> $sampleID
echo "rm $CDT7c" >> $sampleID
echo "rm $CDT7_SCALED" >> $sampleID
echo "rm $PNG7" >> $sampleID
echo "rm $CDT8_SCALED" >> $sampleID
echo "rm $PNG8" >> $sampleID
echo "rm $OUT8a_3000" >> $sampleID
echo "rm $CDT8e_3000" >> $sampleID
echo "rm $CDT8f_3000" >> $sampleID
echo "rm $CDT8_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT23a_3000" >> $sampleID
echo "rm $CDT23e_3000" >> $sampleID
echo "rm $CDT23f_3000" >> $sampleID
echo "rm $CDT23_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT24a_3000" >> $sampleID
echo "rm $CDT24e_3000" >> $sampleID
echo "rm $CDT24f_3000" >> $sampleID
echo "rm $CDT24_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT25a_3000" >> $sampleID
echo "rm $CDT25e_3000" >> $sampleID
echo "rm $CDT25f_3000" >> $sampleID
echo "rm $CDT25_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT26a_3000" >> $sampleID
echo "rm $CDT26_sense_3000_gz" >> $sampleID
echo "rm $CDT26_anti_3000_gz" >> $sampleID
echo "rm $CDT26_anti_3000" >> $sampleID
echo "rm $CDT26_SCALED_COMP_anti_3000" >> $sampleID
echo "rm $OUT32a_3000" >> $sampleID
echo "rm $CDT32_sense_3000_gz" >> $sampleID
echo "rm $CDT32_anti_3000_gz" >> $sampleID
echo "rm $CDT32_sense_3000" >> $sampleID
echo "rm $CDT32_SCALED_COMP_sense_3000" >> $sampleID
echo "rm $OUT28a_3000" >> $sampleID
echo "rm $CDT28e_3000" >> $sampleID
echo "rm $CDT28f_3000" >> $sampleID
echo "rm $CDT28_SCALED_COMP_3000" >> $sampleID
echo "# finish script" >> $sampleID
