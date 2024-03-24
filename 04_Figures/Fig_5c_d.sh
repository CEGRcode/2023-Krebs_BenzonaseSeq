# purpose - make Fig. showing proximal vs. distal hexasome. *this version is based on sort of midpoint of reads with hexasomal size length. ***This version nows uses the FLN-based +1 bedfile, shuffles the bedfile (by row) before during the ratio sort. v3 - showing only top/bottom 1000 sites, but now with MNase ChIP-seq H3K4me3 and matching MNase-seq included. Code (and all new variabes) checked by JEK 240228. This version (v5) no longer includes Benzonase-seq or MNase-seq AND is sorted L vs R based on BX H3K4me3. (v6) - only include anti strand of CoPRO pol II in heatmap.

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
Plus1_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/more_subnucleosomes/fig4c_particles/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet_FullLengthNuc_Only.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig4_5_new_subnucleosomes/H3K4me3_Hex_v6_FLNonly_240318

#set bam library file to BI_rep1 **testing with subsampled master BAM file
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25861_25869_25965_25972_28805_28809_Benz_0sonicCycles_BX_H3K4me3_master.bam
BAM3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/CoPRO/CoPRO_K562_MERGE.bam
BAM4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/19354_19355_NoBenz_10sonicCycles_XO_PolII_master.bam
BAM5=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240215_MNase_data/final_merged_files/SRR6010180_SRR6010175_SRR6010177_SRR7441419_SRR7441420_dedup_MERGE_master.bam

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
#SBATCH --time=6:00:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
module load anaconda
conda activate bioinfo"

#set output file names
BEDFILE_shuffled=$(echo $Plus1_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
BEDFILE_shuffled_a=$(echo $Plus1_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled"}')
Plus1_BEDFILE_2bp=$(echo $BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2bp.bed"}')
Plus1_proximal_hex=$(echo $Plus1_BEDFILE_2bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_proximal_hex.bed"}')
Plus1_distal_hex=$(echo $Plus1_BEDFILE_2bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_distal_hex.bed"}')
Plus1_proximal_hex_a=$(echo $Plus1_proximal_hex | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
Plus1_distal_hex_a=$(echo $Plus1_distal_hex | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$Plus1_proximal_hex_a" | awk -F. '{print $1"_proximal_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$Plus1_proximal_hex_a" | awk -F. '{print $1"_proximal"}')
CDT1_gz=$(echo "$BAM1a""_""$Plus1_proximal_hex_a" | awk -F. '{print $1"_proximal_combined.cdt.gz"}')
CDT1_unzipped=$(echo "$BAM1a""_""$Plus1_proximal_hex_a" | awk -F. '{print $1"_proximal_combined.cdt"}')
OUT2=$(echo "$BAM1a""_""$Plus1_distal_hex_a" | awk -F. '{print $1"_distal_midpoint.out"}')
CDT2=$(echo "$BAM1a""_""$Plus1_distal_hex_a" | awk -F. '{print $1"_distal"}')
CDT2_gz=$(echo "$BAM1a""_""$Plus1_distal_hex_a" | awk -F. '{print $1"_distal_combined.cdt.gz"}')
CDT2_unzipped=$(echo "$BAM1a""_""$Plus1_distal_hex_a" | awk -F. '{print $1"_distal_combined.cdt"}')
target_proximal_hex=$(echo "$BAM1a""_""$Plus1_proximal_hex_a" | awk -F. '{print $1"_RowCount.tab"}')
target_distal_hex=$(echo "$BAM1a""_""$Plus1_distal_hex_a" | awk -F. '{print $1"_RowCount.tab"}')
RATIO_BEDFILE=$(echo "$BEDFILE_shuffled_a" | awk -F. '{print $1"_ratio_proximal_to_distal.bed"}')
RATIO_BEDFILE_200=$(echo $RATIO_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_200bp.bed"}')
RATIO_BEDFILE_2000=$(echo $RATIO_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2000bp.bed"}')
RATIO_BEDFILE_2000_top1000=$(echo $RATIO_BEDFILE_2000| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_top1000.bed"}')
RATIO_BEDFILE_2000_middle1790=$(echo $RATIO_BEDFILE_2000| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_middle1790.bed"}')
RATIO_BEDFILE_2000_bottom1000=$(echo $RATIO_BEDFILE_2000| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_bottom1000.bed"}')
OUT3=$(echo "$BAM1a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_midpoint.out"}')
CDT3=$(echo "$BAM1a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1}')
CDT3b=$(echo "$BAM1a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT3c=$(echo "$BAM1a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_combined.cdt"}')
CDT3_SCALED=$(echo "$BAM1a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled.cdt"}')
SCALE3=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT"}')
SCALE3a=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG3=$(echo "$BAM1a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled.png"}')
SVG3=$(echo "$BAM1a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT3a_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT3d_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite"}')
CDT3e_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT3f_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT3_SCALED_COMP_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT3_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT5=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1.out"}')
CDT5=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1"}')
CDT5_sense_gz=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT5_anti_gz=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT5_sense=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT5_anti=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT5_SCALED_sense=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_sense_scaled.cdt"}')
CDT5_SCALED_anti=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_anti_scaled.cdt"}')
SCALE5=$(echo "$BAM3a" | awk -F. '{print $1"_read1_ForCDT"}')
SCALE5a=$(echo "$BAM3a" | awk -F. '{print $1"_read1_ForCDT_ScalingFactors.out"}')
PNG5_anti=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_scaled_anti.png"}')
SVG5=$(echo "$BAM3a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled_read1_labeled.svg"}')
OUT5a_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite.out"}')
CDT5d_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite"}')
CDT5_sense_2000_gz=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_sense.cdt.gz"}')
CDT5_anti_2000_gz=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt.gz"}')
CDT5_anti_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt"}')
CDT5_SCALED_COMP_anti_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.cdt"}')
SCALED_OUT5_anti_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.tab"}')
BAM4a=$(echo $BAM4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT6=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1.out"}')
CDT6=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1"}')
CDT6_sense_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT6_anti_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT6_sense=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT6_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT6_SCALED_sense=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_sense_scaled.cdt"}')
CDT6_SCALED_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_anti_scaled.cdt"}')
SCALE6=$(echo "$BAM4a" | awk -F. '{print $1"_read1_ForCDT"}')
SCALE6a=$(echo "$BAM4a" | awk -F. '{print $1"_read1_ForCDT_ScalingFactors.out"}')
PNG6_sense=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_scaled_sense.png"}')
PNG6_anti=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_scaled_anti.png"}')
PNG6_merge=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_read1_scaled_merged.png"}')
SVG6=$(echo "$BAM4a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled_read1_labeled.svg"}')
OUT6a_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite.out"}')
CDT6d_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite"}')
CDT6_sense_2000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_sense.cdt.gz"}')
CDT6_anti_2000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt.gz"}')
CDT6_sense_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_sense.cdt"}')
CDT6_anti_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt"}')
CDT6_SCALED_COMP_sense_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_scaled_sense.cdt"}')
CDT6_SCALED_COMP_anti_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.cdt"}')
SCALED_OUT6_sense_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_scaled_sense.tab"}')
SCALED_OUT6_anti_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.tab"}')
SCALED_OUT6_final_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_read1_ForComposite_scaled_final.tab"}')
BAM5a=$(echo $BAM5 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT7=$(echo "$BAM5a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_midpoint.out"}')
CDT7=$(echo "$BAM5a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1}')
CDT7b=$(echo "$BAM5a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT7c=$(echo "$BAM5a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_combined.cdt"}')
CDT7_SCALED=$(echo "$BAM5a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled.cdt"}')
SCALE7=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT"}')
SCALE7a=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG7=$(echo "$BAM5a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled.png"}')
SVG7=$(echo "$BAM5a""_""$RATIO_BEDFILE_200" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT7a_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT7d_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite"}')
CDT7e_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT7f_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT7_SCALED_COMP_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT7_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_top1000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT23a_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT23d_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite"}')
CDT23e_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT23f_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT23_SCALED_COMP_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT23_2000=$(echo "$BAM1a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT25a_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite.out"}')
CDT25d_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite"}')
CDT25_sense_2000_gz=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_sense.cdt.gz"}')
CDT25_anti_2000_gz=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt.gz"}')
CDT25_anti_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt"}')
CDT25_SCALED_COMP_anti_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.cdt"}')
SCALED_OUT25_anti_2000=$(echo "$BAM3a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.tab"}')
OUT26a_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite.out"}')
CDT26d_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite"}')
CDT26_sense_2000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_sense.cdt.gz"}')
CDT26_anti_2000_gz=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt.gz"}')
CDT26_sense_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_sense.cdt"}')
CDT26_anti_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_anti.cdt"}')
CDT26_SCALED_COMP_sense_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_scaled_sense.cdt"}')
CDT26_SCALED_COMP_anti_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.cdt"}')
SCALED_OUT26_sense_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_scaled_sense.tab"}')
SCALED_OUT26_anti_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_scaled_anti.tab"}')
SCALED_OUT26_final_2000=$(echo "$BAM4a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_read1_ForComposite_scaled_final.tab"}')
OUT27a_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT27d_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite"}')
CDT27e_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT27f_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT27_SCALED_COMP_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT27_2000=$(echo "$BAM5a""_""$RATIO_BEDFILE_2000_bottom1000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')

sampleID=H3K4me3_hex_prox_distal_ratio_v6_FLNonly_240318.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#randomize intitial bedfile so that there any no secondary sorts" >> $sampleID
echo "shuf $Plus1_BEDFILE > $BEDFILE_shuffled" >> $sampleID
echo "#expand new bedfile by 2bp so that start and end coodiantes are the midpoint (or dyad) of the nucleosome. Input is the shuffled bedfile." >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2 $BEDFILE_shuffled -o=$Plus1_BEDFILE_2bp" >> $sampleID
echo "#get proximal (upstream) or distal (downstream) 127 bp regions relative to +1 midpoint (dyad)" >> $sampleID
echo "bedtools slop -i $Plus1_BEDFILE_2bp -g $HG19_GENOME -l 72 -r 18 -s > $Plus1_proximal_hex" >> $sampleID
echo "bedtools slop -i $Plus1_BEDFILE_2bp -g $HG19_GENOME -l 18 -r 72 -s > $Plus1_distal_hex" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 --max-insert=127 --min-insert=92 $Plus1_proximal_hex $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 --max-insert=127 --min-insert=92 $Plus1_distal_hex $BAM1" >> $sampleID
echo "#unzip files" >> $sampleID
echo "gunzip -c $CDT1_gz > $CDT1_unzipped" >> $sampleID
echo "gunzip -c $CDT2_gz > $CDT2_unzipped" >> $sampleID
echo "#sum tag in each above CDT file" >> $sampleID
echo "perl $JOB_ROW $CDT1_unzipped $target_proximal_hex" >> $sampleID
echo "perl $JOB_ROW $CDT2_unzipped $target_distal_hex" >> $sampleID
echo "#only keep rows if all IDs match between initial bedfile, prox SN, and distal SN; then remove rows if 0 tags in numerator or denominator, then do division, sort by ratio and print bedfile. There is an exta column (column 7) here." >> $sampleID
echo "paste $BEDFILE_shuffled $target_proximal_hex $target_distal_hex | awk '(\$4==\$8 && \$4==\$10 && \$9!=0 && \$11!=0){print \$0\"\t\"(\$9/\$11)}' | sort -k12,12rn | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $RATIO_BEDFILE" >> $sampleID
echo "#expand new bedfile by 1000 bp and 2000 bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=200 $RATIO_BEDFILE -o=$RATIO_BEDFILE_200" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2000 $RATIO_BEDFILE -o=$RATIO_BEDFILE_2000" >> $sampleID
echo "#take top 2500 sites, middle 5000 sites, and bottom 2500 sites from above above bedfile. Bottom 2500 + ~180 rows that have PolII chunk -> take 2680 but leave last 180" >> $sampleID
echo "cat $RATIO_BEDFILE_2000 | head -1000 > $RATIO_BEDFILE_2000_top1000" >> $sampleID
echo "cat $RATIO_BEDFILE_2000 | tail -1000 > $RATIO_BEDFILE_2000_bottom1000" >> $sampleID
echo "#scale all libraries: options - total tag scaling -t" >> $sampleID
echo "#scale output files: options - total tag scaling -t" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE5 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE6 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE7 $BAM5" >> $sampleID
echo "#do tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 --max-insert=127 --min-insert=92 $RATIO_BEDFILE_200 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT5 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 $RATIO_BEDFILE_200 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT6 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6 $RATIO_BEDFILE_200 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT7 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT7 --max-insert=127 --min-insert=92 $RATIO_BEDFILE_200 $BAM5" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT3b > $OUTPUT/$CDT3c" >> $sampleID
echo "gunzip -c $CDT5_sense_gz > $CDT5_sense" >> $sampleID
echo "gunzip -c $CDT5_anti_gz > $CDT5_anti" >> $sampleID
echo "gunzip -c $CDT6_sense_gz > $CDT6_sense" >> $sampleID
echo "gunzip -c $CDT6_anti_gz > $CDT6_anti" >> $sampleID
echo "gunzip -c $CDT7b > $OUTPUT/$CDT7c" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT3c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_anti --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT5_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_sense --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_anti --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT7_SCALED --scaling-factor=\$(cat $SCALE7a | cut -f2 | tail -1 | awk '{print \$1}') $CDT7c" >> $sampleID
echo "#make heatmap [default black colors used] **change threshold if necessary" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG3 -c ff9900 -p=0.95 $CDT3_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG5_anti -c 0099ff -p=0.95 $CDT5_SCALED_anti" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG6_sense -c 833c0c -p=0.95 $CDT6_SCALED_sense" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG6_anti -c c65911 -p=0.95 $CDT6_SCALED_anti" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG7 -c 000000 -p=0.95 $CDT7_SCALED" >> $sampleID
echo "#merge heatmaps" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation merge-heatmap -o=$PNG6_merge $PNG6_sense $PNG6_anti" >> $sampleID
echo "#label above heatmaps" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG3 --output=$SVG3 --width=2 --font-size=18 --left-label='-0.1' --mid-label=0 --right-label='+0.1' --x-label='Distance from +1 Dyad (kb)' --y-label='8,987 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG5_anti --output=$SVG5 --width=2 --font-size=18 --left-label='-0.1' --mid-label=0 --right-label='+0.1' --x-label='Distance from +1 Dyad (kb)' --y-label='8,987 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG6_merge --output=$SVG6 --width=2 --font-size=18 --left-label='-0.1' --mid-label=0 --right-label='+0.1' --x-label='Distance from +1 Dyad (kb)' --y-label='8,987 Nucleosome positions sorted by ratio'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG7 --output=$SVG7 --width=2 --font-size=18 --left-label='-0.1' --mid-label=0 --right-label='+0.1' --x-label='Distance from +1 Dyad (kb)' --y-label='8,987 Nucleosome positions sorted by ratio'" >> $sampleID
echo "#prep files for comp plots at top 1000 sites of bedfile" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3d_2000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3a_2000 --max-insert=127 --min-insert=92 $RATIO_BEDFILE_2000_top1000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT5d_2000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5a_2000 $RATIO_BEDFILE_2000_top1000 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT6d_2000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6a_2000 $RATIO_BEDFILE_2000_top1000 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT7d_2000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT7a_2000 --max-insert=127 --min-insert=92 $RATIO_BEDFILE_2000_top1000 $BAM5" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT3e_2000 > $OUTPUT/$CDT3f_2000" >> $sampleID
echo "gunzip -c $CDT5_anti_2000_gz > $CDT5_anti_2000" >> $sampleID
echo "gunzip -c $CDT6_sense_2000_gz > $CDT6_sense_2000" >> $sampleID
echo "gunzip -c $CDT6_anti_2000_gz > $CDT6_anti_2000" >> $sampleID
echo "gunzip -c $CDT7e_2000 > $OUTPUT/$CDT7f_2000" >> $sampleID
echo "#scale CDT files data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_2000 --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT3f_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_COMP_anti_2000 --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT5_anti_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_COMP_sense_2000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_sense_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_COMP_anti_2000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_anti_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT7_SCALED_COMP_2000 --scaling-factor=\$(cat $SCALE7a | cut -f2 | tail -1 | awk '{print \$1}') $CDT7f_2000" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_2000 $SCALED_OUT3_2000" >> $sampleID
echo "perl $JOB $CDT5_SCALED_COMP_anti_2000 $SCALED_OUT5_anti_2000" >> $sampleID
echo "perl $JOB $CDT6_SCALED_COMP_sense_2000 $SCALED_OUT6_sense_2000" >> $sampleID
echo "perl $JOB $CDT6_SCALED_COMP_anti_2000 $SCALED_OUT6_anti_2000" >> $sampleID
echo "perl $JOB $CDT7_SCALED_COMP_2000 $SCALED_OUT7_2000" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library. FILE NOT made for subfigs 3 and 4 has has 2 stranded are COMBINDED." >> $sampleID
echo "cat $SCALED_OUT6_sense_2000 $SCALED_OUT6_anti_2000 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT6_final_2000" >> $sampleID
echo "#prep files for comp plots at bottom 1000 sites of bedfile" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT23d_2000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT23a_2000 --max-insert=127 --min-insert=92 $RATIO_BEDFILE_2000_bottom1000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT25d_2000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT25a_2000 $RATIO_BEDFILE_2000_bottom1000 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT26d_2000 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT26a_2000 $RATIO_BEDFILE_2000_bottom1000 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT27d_2000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT27a_2000 --max-insert=127 --min-insert=92 $RATIO_BEDFILE_2000_bottom1000 $BAM5" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT23e_2000 > $OUTPUT/$CDT23f_2000" >> $sampleID
echo "gunzip -c $CDT25_anti_2000_gz > $CDT25_anti_2000" >> $sampleID
echo "gunzip -c $CDT26_sense_2000_gz > $CDT26_sense_2000" >> $sampleID
echo "gunzip -c $CDT26_anti_2000_gz > $CDT26_anti_2000" >> $sampleID
echo "gunzip -c $CDT27e_2000 > $OUTPUT/$CDT27f_2000" >> $sampleID
echo "#scale CDT files data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT23_SCALED_COMP_2000 --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT23f_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT25_SCALED_COMP_anti_2000 --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT25_anti_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT26_SCALED_COMP_sense_2000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT26_sense_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT26_SCALED_COMP_anti_2000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT26_anti_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT27_SCALED_COMP_2000 --scaling-factor=\$(cat $SCALE7a | cut -f2 | tail -1 | awk '{print \$1}') $CDT27f_2000" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT23_SCALED_COMP_2000 $SCALED_OUT23_2000" >> $sampleID
echo "perl $JOB $CDT25_SCALED_COMP_anti_2000 $SCALED_OUT25_anti_2000" >> $sampleID
echo "perl $JOB $CDT26_SCALED_COMP_sense_2000 $SCALED_OUT26_sense_2000" >> $sampleID
echo "perl $JOB $CDT26_SCALED_COMP_anti_2000 $SCALED_OUT26_anti_2000" >> $sampleID
echo "perl $JOB $CDT27_SCALED_COMP_2000 $SCALED_OUT27_2000" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library. FILE NOT made for subfigs 3 and 4 has has 2 stranded are COMBINDED." >> $sampleID
echo "cat $SCALED_OUT26_sense_2000 $SCALED_OUT26_anti_2000 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT26_final_2000" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1_gz" >> $sampleID
echo "rm $CDT1_unzipped" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2_gz" >> $sampleID
echo "rm $CDT2_unzipped" >> $sampleID
echo "rm $target_proximal_hex" >> $sampleID
echo "rm $target_distal_hex" >> $sampleID
echo "rm $RATIO_BEDFILE_2000" >> $sampleID
echo "rm $RATIO_BEDFILE_2000" >> $sampleID
echo "rm $OUT3" >> $sampleID
echo "rm $CDT3b" >> $sampleID
echo "rm $CDT3c" >> $sampleID
echo "rm $CDT3_SCALED" >> $sampleID
echo "rm $SCALE3a" >> $sampleID
echo "rm $PNG3" >> $sampleID
echo "rm $OUT3a_2000" >> $sampleID
echo "rm $CDT3e_2000" >> $sampleID
echo "rm $CDT3f_2000" >> $sampleID
echo "rm $CDT3_SCALED_COMP_2000" >> $sampleID
echo "rm $OUT5" >> $sampleID
echo "rm $CDT5_sense_gz" >> $sampleID
echo "rm $CDT5_anti_gz" >> $sampleID
echo "rm $CDT5_sense" >> $sampleID
echo "rm $CDT5_anti" >> $sampleID
echo "rm $CDT5_SCALED_sense" >> $sampleID
echo "rm $CDT5_SCALED_anti" >> $sampleID
echo "rm $SCALE5a" >> $sampleID
echo "rm $PNG5_anti" >> $sampleID
echo "rm $OUT5a_2000" >> $sampleID
echo "rm $CDT5_sense_2000_gz" >> $sampleID
echo "rm $CDT5_anti_2000_gz" >> $sampleID
echo "rm $CDT5_anti_2000" >> $sampleID
echo "rm $CDT5_SCALED_COMP_anti_2000" >> $sampleID
echo "rm $CDT6_sense_gz" >> $sampleID
echo "rm $CDT6_anti_gz" >> $sampleID
echo "rm $CDT6_sense" >> $sampleID
echo "rm $CDT6_anti" >> $sampleID
echo "rm $CDT6_SCALED_sense" >> $sampleID
echo "rm $CDT6_SCALED_anti" >> $sampleID
echo "rm $SCALE6a" >> $sampleID
echo "rm $PNG6_sense" >> $sampleID
echo "rm $PNG6_anti" >> $sampleID
echo "rm $PNG6_merge" >> $sampleID
echo "rm $OUT6a_2000" >> $sampleID
echo "rm $CDT6_sense_2000_gz" >> $sampleID
echo "rm $CDT6_anti_2000_gz" >> $sampleID
echo "rm $CDT6_sense_2000" >> $sampleID
echo "rm $CDT6_anti_2000" >> $sampleID
echo "rm $CDT6_SCALED_COMP_sense_2000" >> $sampleID
echo "rm $CDT6_SCALED_COMP_anti_2000" >> $sampleID
echo "rm $SCALED_OUT6_sense_2000" >> $sampleID
echo "rm $SCALED_OUT6_anti_2000" >> $sampleID
echo "rm $OUT7" >> $sampleID
echo "rm $CDT7b" >> $sampleID
echo "rm $CDT7c" >> $sampleID
echo "rm $CDT7_SCALED" >> $sampleID
echo "rm $SCALE7a" >> $sampleID
echo "rm $PNG7" >> $sampleID
echo "rm $OUT7a_2000" >> $sampleID
echo "rm $CDT7e_2000" >> $sampleID
echo "rm $CDT7f_2000" >> $sampleID
echo "rm $CDT7_SCALED_COMP_2000" >> $sampleID
echo "rm $OUT23a_2000" >> $sampleID
echo "rm $CDT23e_2000" >> $sampleID
echo "rm $CDT23f_2000" >> $sampleID
echo "rm $CDT23_SCALED_COMP_2000" >> $sampleID
echo "rm $OUT25a_2000" >> $sampleID
echo "rm $CDT25_sense_2000_gz" >> $sampleID
echo "rm $CDT25_anti_2000_gz" >> $sampleID
echo "rm $CDT25_anti_2000" >> $sampleID
echo "rm $CDT25_SCALED_COMP_anti_2000" >> $sampleID
echo "rm $OUT26a_2000" >> $sampleID
echo "rm $CDT26_sense_2000_gz" >> $sampleID
echo "rm $CDT26_anti_2000_gz" >> $sampleID
echo "rm $CDT26_sense_2000" >> $sampleID
echo "rm $CDT26_anti_2000" >> $sampleID
echo "rm $CDT26_SCALED_COMP_sense_2000" >> $sampleID
echo "rm $CDT26_SCALED_COMP_anti_2000" >> $sampleID
echo "rm $SCALED_OUT26_sense_2000" >> $sampleID
echo "rm $SCALED_OUT26_anti_2000" >> $sampleID
echo "rm $OUT27a_2000" >> $sampleID
echo "rm $CDT27e_2000" >> $sampleID
echo "rm $CDT27f_2000" >> $sampleID
echo "rm $CDT27_SCALED_COMP_2000" >> $sampleID
echo "# finish script" >> $sampleID
