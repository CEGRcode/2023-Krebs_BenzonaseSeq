# purpose - calculate density for all +1 nucleosomes (11,714 nucleosomes) with all active histone mod(s) and their respective histone. Code checked 240117. (v7) - edges of bedfiles changed so that +/-73 bp from dyad of +1 nucleosome.
# usage
# qq
#
# example
#
# 'qq'

#set bedfiles with  Plus1 (downstream) nucleosome positions based on non-redundant bedfile-based calls
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/nre_fig3_240219/fig3d_density_128to164bp_v7_output_240320

#set bam files
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25861_25869_25965_25972_28805_28809_Benz_0sonicCycles_BX_H3K4me3_master.bam
BAM2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25862_25870_25966_25973_Benz_0sonicCycles_BX_H3K9Ac_master.bam
BAM3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25858_25866_25962_25969_28806_28810_Benz_0sonicCycles_BX_H3K27Ac_master.bam
BAM4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/25860_25868_25964_25971_28804_28808_Benz_0sonicCycles_BX_H3_master.bam
BAM5=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/28453_28794_28461_28799_Benz_0sonicCycles_BX_H2A_Z_master.bam
BAM6=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_ChIPs/MERGED_datasets/28454_28795_28462_28800_Benz_0sonicCycles_BX_H2B_master.bam

#set genome and human.hg19.genome file
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/GENOMES/hg19.fa
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl
PLOT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig5_density/bedfile_RNAvalues_modified_Seaborn_240303/violin_plots_mod2_240303.py

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
#SBATCH --time=1:00:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
conda activate plot"

#set output file names
BEDFILE_2bp=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2bp.bed"}')
Plus1_proximal_half=$(echo $BEDFILE_2bp| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_proximal_half.bed"}')
Plus1_distal_half=$(echo $BEDFILE_2bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_distal_half.bed"}')
Plus1_proximal_half_a=$(echo $Plus1_proximal_half | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
Plus1_distal_half_a=$(echo $Plus1_distal_half | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1b=$(echo $BAM1 | rev | cut -d"X" -f1 | rev | awk '{print substr($1,2);}' | awk -F"." '{print $1}')
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM2b=$(echo $BAM2 | rev | cut -d"X" -f1 | rev | awk '{print substr($1,2);}' | awk -F"." '{print $1}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM3b=$(echo $BAM3 | rev | cut -d"X" -f1 | rev | awk '{print substr($1,2);}' | awk -F"." '{print $1}')
BAM4a=$(echo $BAM4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM4b=$(echo $BAM4 | rev | cut -d"X" -f1 | rev | awk '{print substr($1,2);}' | awk -F"." '{print $1}')
BAM5a=$(echo $BAM5 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM5b=$(echo $BAM5 | rev | cut -d"X" -f1 | rev | awk '{print substr($1,2);}' | awk -F"." '{print $1}')
BAM6a=$(echo $BAM6 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM6b=$(echo $BAM6 | rev | cut -d"X" -f1 | rev | awk '{print substr($1,2);}' | awk -F"." '{print $1}')
OUT1=$(echo "$BAM1a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1.out"}')
CDT1=$(echo "$BAM1a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1"}')
OUT2=$(echo "$BAM1a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1.out"}')
CDT2=$(echo "$BAM1a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1"}')
OUT3=$(echo "$BAM2a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1.out"}')
CDT3=$(echo "$BAM2a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1"}')
OUT4=$(echo "$BAM2a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1.out"}')
CDT4=$(echo "$BAM2a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1"}')
OUT5=$(echo "$BAM3a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1.out"}')
CDT5=$(echo "$BAM3a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1"}')
OUT6=$(echo "$BAM3a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1.out"}')
CDT6=$(echo "$BAM3a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1"}')
OUT7=$(echo "$BAM4a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1.out"}')
CDT7=$(echo "$BAM4a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1"}')
OUT8=$(echo "$BAM4a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1.out"}')
CDT8=$(echo "$BAM4a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1"}')
OUT9=$(echo "$BAM5a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1.out"}')
CDT9=$(echo "$BAM5a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1"}')
OUT10=$(echo "$BAM5a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1.out"}')
CDT10=$(echo "$BAM5a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1"}')
OUT11=$(echo "$BAM6a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1.out"}')
CDT11=$(echo "$BAM6a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1"}')
OUT12=$(echo "$BAM6a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1.out"}')
CDT12=$(echo "$BAM6a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1"}')
CDT1_sense_gz=$(echo "$BAM1a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT3_sense_gz=$(echo "$BAM2a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT5_sense_gz=$(echo "$BAM3a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT7_sense_gz=$(echo "$BAM4a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT9_sense_gz=$(echo "$BAM5a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT11_sense_gz=$(echo "$BAM6a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT2_anti_gz=$(echo "$BAM1a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT4_anti_gz=$(echo "$BAM2a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT6_anti_gz=$(echo "$BAM3a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT8_anti_gz=$(echo "$BAM4a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT10_anti_gz=$(echo "$BAM5a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT12_anti_gz=$(echo "$BAM6a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT1_sense=$(echo "$BAM1a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT3_sense=$(echo "$BAM2a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT5_sense=$(echo "$BAM3a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT7_sense=$(echo "$BAM4a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT9_sense=$(echo "$BAM5a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT11_sense=$(echo "$BAM6a""_""$Plus1_proximal_half" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT2_anti=$(echo "$BAM1a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT4_anti=$(echo "$BAM2a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT6_anti=$(echo "$BAM3a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT8_anti=$(echo "$BAM4a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT10_anti=$(echo "$BAM5a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT12_anti=$(echo "$BAM6a""_""$Plus1_distal_half" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT1_sense_sum=$(echo "$BAM1a""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_sum.tsv"}')
CDT3_sense_sum=$(echo "$BAM2a""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_sum.tsv"}')
CDT5_sense_sum=$(echo "$BAM3a""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_sum.tsv"}')
CDT7_sense_sum=$(echo "$BAM4a""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_sum.tsv"}')
CDT9_sense_sum=$(echo "$BAM5a""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_sum.tsv"}')
CDT11_sense_sum=$(echo "$BAM6a""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_sum.tsv"}')
CDT2_anti_sum=$(echo "$BAM1a""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_sum.tsv"}')
CDT4_anti_sum=$(echo "$BAM2a""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_sum.tsv"}')
CDT6_anti_sum=$(echo "$BAM3a""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_sum.tsv"}')
CDT8_anti_sum=$(echo "$BAM4a""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_sum.tsv"}')
CDT10_anti_sum=$(echo "$BAM5a""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_sum.tsv"}')
CDT12_anti_sum=$(echo "$BAM6a""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_sum.tsv"}')
CDT1_CDT7_sense_bins_sum=$(echo "$BAM1b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_bins_sum.tsv"}')
CDT2_CDT8_anti_bins_sum=$(echo "$BAM1b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_bins_sum.tsv"}')
CDT3_CDT7_sense_bins_sum=$(echo "$BAM2b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_bins_sum.tsv"}')
CDT4_CDT8_anti_bins_sum=$(echo "$BAM2b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_bins_sum.tsv"}')
CDT5_CDT7_sense_bins_sum=$(echo "$BAM3b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_bins_sum.tsv"}')
CDT6_CDT8_anti_bins_sum=$(echo "$BAM3b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_bins_sum.tsv"}')
CDT9_CDT11_sense_bins_sum=$(echo "$BAM5b""_""$BAM6b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_bins_sum.tsv"}')
CDT10_CDT12_anti_bins_sum=$(echo "$BAM5b""_""$BAM6b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_bins_sum.tsv"}')
CDT1_CDT7_sense_density=$(echo "$BAM1b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density.tsv"}')
CDT2_CDT8_anti_density=$(echo "$BAM1b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density.tsv"}')
CDT3_CDT7_sense_density=$(echo "$BAM2b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density.tsv"}')
CDT4_CDT8_anti_density=$(echo "$BAM2b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density.tsv"}')
CDT5_CDT7_sense_density=$(echo "$BAM3b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density.tsv"}')
CDT6_CDT8_anti_density=$(echo "$BAM3b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density.tsv"}')
CDT9_CDT11_sense_density=$(echo "$BAM5b""_""$BAM6b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density.tsv"}')
CDT10_CDT12_anti_density=$(echo "$BAM5b""_""$BAM6b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density.tsv"}')
CDT1_CDT7_sense_density_ID=$(echo "$BAM1b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density_ID.tsv"}')
CDT2_CDT8_anti_density_ID=$(echo "$BAM1b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density_ID.tsv"}')
CDT3_CDT7_sense_density_ID=$(echo "$BAM2b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density_ID.tsv"}')
CDT4_CDT8_anti_density_ID=$(echo "$BAM2b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density_ID.tsv"}')
CDT5_CDT7_sense_density_ID=$(echo "$BAM3b""_""$BAM4b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density_ID.tsv"}')
CDT6_CDT8_anti_density_ID=$(echo "$BAM3b""_""$BAM4b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density_ID.tsv"}')
CDT9_CDT11_sense_density_ID=$(echo "$BAM5b""_""$BAM6b""_""$Plus1_proximal_half" | awk -F. '{print $1"_sense_density_ID.tsv"}')
CDT10_CDT12_anti_density_ID=$(echo "$BAM5b""_""$BAM6b""_""$Plus1_distal_half" | awk -F. '{print $1"_anti_density_ID.tsv"}')
final_density_input=$(echo "histoneMods" | awk -F. '{print $1"_density_file.tsv"}')
density_file=$(echo "histoneMods" | awk -F. '{print $1"_density_file_for_correlation.tsv"}')
SVG=$(echo "histoneMods" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1".svg"}')

sampleID=Fig3d_HistoneMods_density_v7_128to164bp_240320.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#expand bedfiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2 $BEDFILE -o=$BEDFILE_2bp" >> $sampleID
echo "#get proximal (upstream) or distal (downstream) half of nucleosome relative to +1 midpoint (dyad)" >> $sampleID
echo "bedtools slop -i $BEDFILE_2bp -g $HG19_GENOME -l 72 -r -2 -s > $Plus1_proximal_half" >> $sampleID
echo "bedtools slop -i $BEDFILE_2bp -g $HG19_GENOME -l -2 -r 72 -s > $Plus1_distal_half" >> $sampleID
echo "#run tag pileup" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT1 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 --min-insert=128 --max-insert=164 $Plus1_proximal_half $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT2 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 --min-insert=128 --max-insert=164 $Plus1_distal_half $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT3 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 --min-insert=128 --max-insert=164 $Plus1_proximal_half $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT4 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 --min-insert=128 --max-insert=164 $Plus1_distal_half $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT5 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 --min-insert=128 --max-insert=164 $Plus1_proximal_half $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT6 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6 --min-insert=128 --max-insert=164 $Plus1_distal_half $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT7 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT7 --min-insert=128 --max-insert=164 $Plus1_proximal_half $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT8 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT8 --min-insert=128 --max-insert=164 $Plus1_distal_half $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT9 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT9 --min-insert=128 --max-insert=164 $Plus1_proximal_half $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT10 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT10 --min-insert=128 --max-insert=164 $Plus1_distal_half $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT11 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT11 --min-insert=128 --max-insert=164 $Plus1_proximal_half $BAM6" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT12 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT12 --min-insert=128 --max-insert=164 $Plus1_distal_half $BAM6" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1_sense_gz > $CDT1_sense" >> $sampleID
echo "gunzip -c $CDT2_anti_gz > $CDT2_anti" >> $sampleID
echo "gunzip -c $CDT3_sense_gz > $CDT3_sense" >> $sampleID
echo "gunzip -c $CDT4_anti_gz > $CDT4_anti" >> $sampleID
echo "gunzip -c $CDT5_sense_gz > $CDT5_sense" >> $sampleID
echo "gunzip -c $CDT6_anti_gz > $CDT6_anti" >> $sampleID
echo "gunzip -c $CDT7_sense_gz > $CDT7_sense" >> $sampleID
echo "gunzip -c $CDT8_anti_gz > $CDT8_anti" >> $sampleID
echo "gunzip -c $CDT9_sense_gz > $CDT9_sense" >> $sampleID
echo "gunzip -c $CDT10_anti_gz > $CDT10_anti" >> $sampleID
echo "gunzip -c $CDT11_sense_gz > $CDT11_sense" >> $sampleID
echo "gunzip -c $CDT12_anti_gz > $CDT12_anti" >> $sampleID
echo "#sum the number of tags by each row" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT1_sense_sum -r=1 $CDT1_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT2_anti_sum -r=1 $CDT2_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT3_sense_sum -r=1 $CDT3_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT4_anti_sum -r=1 $CDT4_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT5_sense_sum -r=1 $CDT5_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT6_anti_sum -r=1 $CDT6_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT7_sense_sum -r=1 $CDT7_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT8_anti_sum -r=1 $CDT8_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT9_sense_sum -r=1 $CDT9_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT10_anti_sum -r=1 $CDT10_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT11_sense_sum -r=1 $CDT11_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT12_anti_sum -r=1 $CDT12_anti" >> $sampleID
echo "#join CDT1 and CDT2 sense file (or anti files) and remove header line" >> $sampleID
echo "paste $CDT1_sense_sum $CDT7_sense_sum | awk 'NR>1' > $CDT1_CDT7_sense_bins_sum" >> $sampleID
echo "paste $CDT2_anti_sum $CDT8_anti_sum | awk 'NR>1' > $CDT2_CDT8_anti_bins_sum" >> $sampleID
echo "paste $CDT3_sense_sum $CDT7_sense_sum | awk 'NR>1' > $CDT3_CDT7_sense_bins_sum" >> $sampleID
echo "paste $CDT4_anti_sum $CDT8_anti_sum | awk 'NR>1' > $CDT4_CDT8_anti_bins_sum" >> $sampleID
echo "paste $CDT5_sense_sum $CDT7_sense_sum | awk 'NR>1' > $CDT5_CDT7_sense_bins_sum" >> $sampleID
echo "paste $CDT6_anti_sum $CDT8_anti_sum | awk 'NR>1' > $CDT6_CDT8_anti_bins_sum" >> $sampleID
echo "paste $CDT9_sense_sum $CDT11_sense_sum | awk 'NR>1' > $CDT9_CDT11_sense_bins_sum" >> $sampleID
echo "paste $CDT10_anti_sum $CDT12_anti_sum | awk 'NR>1' > $CDT10_CDT12_anti_bins_sum" >> $sampleID
echo "#calculate density" >> $sampleID
echo "cat $CDT1_CDT7_sense_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT1_CDT7_sense_density" >> $sampleID
echo "cat $CDT2_CDT8_anti_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT2_CDT8_anti_density" >> $sampleID
echo "cat $CDT3_CDT7_sense_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT3_CDT7_sense_density" >> $sampleID
echo "cat $CDT4_CDT8_anti_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT4_CDT8_anti_density" >> $sampleID
echo "cat $CDT5_CDT7_sense_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT5_CDT7_sense_density" >> $sampleID
echo "cat $CDT6_CDT8_anti_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT6_CDT8_anti_density" >> $sampleID
echo "cat $CDT9_CDT11_sense_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT9_CDT11_sense_density" >> $sampleID
echo "cat $CDT10_CDT12_anti_bins_sum | awk 'BEGIN {OFS=\"\t\"}{x = \$2; y = \$4; z = log((x+1)/(y+1))log(2); printf \"%.3f\n\", z}' > $CDT10_CDT12_anti_density" >> $sampleID
echo "#add column with ID of file" >> $sampleID
echo "cat $CDT1_CDT7_sense_density | awk '{print \$1\"\t\"\"$BAM1b\"\"/\"\"$BAM4b\"\"_proximal_sense\"}' > $CDT1_CDT7_sense_density_ID" >> $sampleID
echo "cat $CDT2_CDT8_anti_density | awk '{print \$1\"\t\"\"$BAM1b\"\"/\"\"$BAM4b\"\"_distal_anti\"}' > $CDT2_CDT8_anti_density_ID" >> $sampleID
echo "cat $CDT3_CDT7_sense_density | awk '{print \$1\"\t\"\"$BAM2b\"\"/\"\"$BAM4b\"\"_proximal_sense\"}' > $CDT3_CDT7_sense_density_ID" >> $sampleID
echo "cat $CDT4_CDT8_anti_density | awk '{print \$1\"\t\"\"$BAM2b\"\"/\"\"$BAM4b\"\"_distal_anti\"}' > $CDT4_CDT8_anti_density_ID" >> $sampleID
echo "cat $CDT5_CDT7_sense_density | awk '{print \$1\"\t\"\"$BAM3b\"\"/\"\"$BAM4b\"\"_proximal_sense\"}' > $CDT5_CDT7_sense_density_ID" >> $sampleID
echo "cat $CDT6_CDT8_anti_density | awk '{print \$1\"\t\"\"$BAM3b\"\"/\"\"$BAM4b\"\"_distal_anti\"}' > $CDT6_CDT8_anti_density_ID" >> $sampleID
echo "cat $CDT9_CDT11_sense_density | awk '{print \$1\"\t\"\"$BAM5b\"\"/\"\"$BAM6b\"\"_proximal_sense\"}' > $CDT9_CDT11_sense_density_ID" >> $sampleID
echo "cat $CDT10_CDT12_anti_density | awk '{print \$1\"\t\"\"$BAM5b\"\"/\"\"$BAM6b\"\"_distal_anti\"}' > $CDT10_CDT12_anti_density_ID" >> $sampleID
echo "#make file of all data" >> $sampleID
echo "cat $CDT1_CDT7_sense_density_ID $CDT2_CDT8_anti_density_ID $CDT3_CDT7_sense_density_ID $CDT4_CDT8_anti_density_ID $CDT5_CDT7_sense_density_ID $CDT6_CDT8_anti_density_ID $CDT9_CDT11_sense_density_ID $CDT10_CDT12_anti_density_ID > $final_density_input" >> $sampleID
echo "#make violin plot wity modified violin_plots_mod.py" >> $sampleID
echo "python $PLOT -i $final_density_input -o $SVG" >> $sampleID
echo "#make file to check correaltion of all densities" >> $sampleID
echo "paste $BEDFILE $CDT1_CDT7_sense_density $CDT2_CDT8_anti_density $CDT3_CDT7_sense_density $CDT4_CDT8_anti_density $CDT5_CDT7_sense_density $CDT6_CDT8_anti_density $CDT9_CDT11_sense_density $CDT10_CDT12_anti_density > $density_file" >> $sampleID
echo "# finish script" >> $sampleID
