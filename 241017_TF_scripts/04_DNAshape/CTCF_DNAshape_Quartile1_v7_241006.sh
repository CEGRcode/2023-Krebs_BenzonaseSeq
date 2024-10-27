# purpose - iake bedfiles of all quartiles from TF/nuc ratio script -> make bedfiles of all quartiles -> runs DNA shape. v5 - 240922; v6 - more code

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
QUARTILE1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/CTCF_NucOccupancy_settings_pipeline_MA1929_1_240910/MA1929_1_final_1000bp_intersected_164bp_category1_1000bp.bed
MEME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240830_TFBS/MEME_files_240904/MA1929.1.meme

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/03_DNAshape_240922/CTCF_DNAshape_MA1929_1_Quartile1_v7_241006

#set scriptmanager and other jobs
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
SMOOTH3=TF_pipeline_jobs/smoothing_240813.py
EXTRACT=TF_pipeline_jobs/extract_row_number_240817.py
MASKED=DNAshape_jobs/masked_region_DNAshape_241006.py
MAX_MIN_SCALE=DNAshape_jobs/max_min_scale_v2_241006.py
FORMAT=DNAshape_jobs/format_240922.py
FINAL=DNAshape_jobs/final_240922.py

#set genome and human.hg38.genome file
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/hg38_genome/hg38.fa
HG38_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240830_TFBS/FILES/ENCFF667IGK.tsv

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

module load anaconda #configures shell to use conda activate
conda activate plot"

#set output file names
DNAshape=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HelT=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_HelT.out"}')
MGW=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_MGW.out"}')
PropT=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_PropT.out"}')
Roll=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_Roll.out"}')
HelT_3=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_HelT_smoothed_3bp.tab"}')
MGW_3=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MGW_smoothed_3bp.tab"}')
PropT_3=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_PropT_smoothed_3bp.tab"}')
Roll_3=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_Roll_smoothed_3bp.tab"}')
HelT_scale=$(echo "$HelT_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
MGW_scale=$(echo "$MGW_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
PropT_scale=$(echo "$PropT_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
Roll_scale=$(echo "$Roll_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
HelT_final=$(echo "01_""$HelT" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
MGW_final=$(echo "02_""$MGW" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
PropT_final=$(echo "03_""$PropT" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
Roll_final=$(echo "04_""$Roll" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
FILE_final=$(echo "$QUARTILE1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final_file.tab"}')
NT_count=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NT_count.tab"}')
MASKED_region=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_masked.tab"}')

sampleID=CTCF_DNAshape_Quartile1_v7_241006.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#determine DNA shape with scriptmanager" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed --all --avg-composite -o=$DNAshape $GENOME $QUARTILE1" >> $sampleID
echo "#apply 3 bp smoothing" >> $sampleID
echo "python3 $SMOOTH3 $HelT $HelT_3" >> $sampleID
echo "python3 $SMOOTH3 $MGW $MGW_3" >> $sampleID
echo "python3 $SMOOTH3 $PropT $PropT_3" >> $sampleID
echo "python3 $SMOOTH3 $Roll $Roll_3" >> $sampleID
echo "#extract number of NTs from MEME file" >> $sampleID
echo "python3 $EXTRACT $MEME $NT_count" >> $sampleID
echo "#determine the 5' and 3' boundaries of the motif masked region relative to the center column of tab files at column 256" >> $sampleID
echo "python3 $MASKED $NT_count $MASKED_region" >> $sampleID
echo "#determine max scale for +/- 2bp arbitary units" >> $sampleID
echo "python3 $MAX_MIN_SCALE $HelT_3 $MASKED_region $HelT_scale" >> $sampleID
echo "python3 $MAX_MIN_SCALE $MGW_3 $MASKED_region $MGW_scale" >> $sampleID
echo "python3 $MAX_MIN_SCALE $PropT_3 $MASKED_region $PropT_scale" >> $sampleID
echo "python3 $MAX_MIN_SCALE $Roll_3 $MASKED_region $Roll_scale" >> $sampleID
echo "#change file name of OUT file and swap signs (if most of values are negative)" >> $sampleID
echo "python3 $FORMAT $HelT_scale $HelT $HelT_final" >> $sampleID
echo "python3 $FORMAT $MGW_scale $MGW $MGW_final" >> $sampleID
echo "python3 $FORMAT $PropT_scale $PropT $PropT_final" >> $sampleID
echo "python3 $FORMAT $Roll_scale $Roll $Roll_final" >> $sampleID
echo "#make file of scale and indication if values were swapped" >> $sampleID
echo "python3 $FINAL $HelT_scale $MGW_scale $PropT_scale $Roll_scale $FILE_final" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $HelT" >> $sampleID
echo "rm $MGW" >> $sampleID
echo "rm $PropT" >> $sampleID
echo "rm $Roll" >> $sampleID
echo "rm $HelT_3" >> $sampleID
echo "rm $MGW_3" >> $sampleID
echo "rm $PropT_3" >> $sampleID
echo "rm $Roll_3" >> $sampleID
echo "#script DONE" >> $sampleID
