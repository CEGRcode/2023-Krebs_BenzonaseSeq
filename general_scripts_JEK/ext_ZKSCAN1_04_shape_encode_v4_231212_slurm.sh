# purpose - intersect bedfile of all sites for a motif with encode-called peaks (Bed narrowPeak bedfile). Then sort by column 7 of the encode peak. Take the bottom 1/2 of 'bound sites' and look at BNase-seq at the 'lowly bound sites'

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

#set bedfiles
LOWLY_BOUND=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/extended_ZKSCAN1_motif_231211/benzonase_output/MA1585_1_final_1000bp_intersected_lowlyBound_1000bp_motif1_meme_1000bp_shifted_final_1000bp_intersected_lowlyBound_1000bp.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/extended_ZKSCAN1_motif_231211/04_DNAshape_output

#set bam library file to BI_rep1 **testing with subsampled master BAM file
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/K562_benzonase-seq_master.bam

#set blacklist and .genome file
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/GENOMES/hg19.fa
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
#SBATCH --time=2:00:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

#set output file names
DNAshape_lowlyBound=$(echo $LOWLY_BOUND | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')

sampleID=ext_ZKSCAN1_DNAshape_lowlyBound_v4_231212.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#determine DNA shape with scriptmanager" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed --all --avg-composite -o=$DNAshape_lowlyBound $GENOME $LOWLY_BOUND" >> $sampleID
echo "#script DONE" >> $sampleID
