# purpose - calculate GRO-cap 
# usage
# qq
#
# example
#
# 'qq'

#set bedfiles with  Plus1 (downstream) nucleosome positions based on non-redundant bedfile-based calls
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240321_GROcap_RNA/GRO_cap_RNA_tota_v2_output_240321

#set bam files
BAM1=/storage/group/bfp2/default/hxc585_HainingChen/2023_Chen_PIC3/data/RNA_seq/K562_Gro-cap_ENCFF028THC.bam

#set genome and human.hg19.genome file
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/GENOMES/hg19.fa
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl
JOB_ROW=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig6_Subnucleosomes/job/sum_Row_CDT.pl
PLOT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240314_pausing_index/violin_plots_mod240315.py

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
#SBATCH --time=0:30:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
conda activate plot"

#set output file names
BEDFILE_400bp=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_400bp.bed"}')
BEDFILE_400bp_a=$(echo $BEDFILE_400bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$BEDFILE_400bp" | awk -F. '{print $1"_read1.out"}')
CDT1=$(echo "$BAM1a""_""$BEDFILE_400bp" | awk -F. '{print $1"_read1"}')
CDT1_sense_gz=$(echo "$BAM1a""_""$BEDFILE_400bp" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT1_anti_gz=$(echo "$BAM1a""_""$BEDFILE_400bp" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT1_sense=$(echo "$BAM1a""_""$BEDFILE_400bp" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT1_sense_sum=$(echo "$BAM1a""_""$BEDFILE_400bp" | awk -F. '{print $1"_sense_sum.tsv"}')
CDT1_sense_bins_sum=$(echo "$BAM1a" | awk -F. '{print $1"_sense_bins_sum.tsv"}')
RNAtotal=$(echo "$BEDFILE" | awk -F. '{print $1"_GROcap_total.tsv"}')

sampleID=GROcap_RNA_total_v2.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#expand bedfiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=400 $BEDFILE -o=$BEDFILE_400bp" >> $sampleID
echo "#run tag pileup" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -1 -5 -z --output-matrix=$CDT1 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $BEDFILE_400bp $BAM1" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1_sense_gz > $CDT1_sense" >> $sampleID
echo "#sum the number of tags by each row" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT1_sense_sum -r=1 $CDT1_sense" >> $sampleID
echo "#join CDT1 and CDT2 anti files and remove header line" >> $sampleID
echo "paste $CDT1_sense_sum  | awk 'NR>1' > $CDT1_sense_bins_sum" >> $sampleID
