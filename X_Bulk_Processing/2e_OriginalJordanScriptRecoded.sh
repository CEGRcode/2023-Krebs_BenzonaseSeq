#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24GB
#SBATCH --time=2:00:00
#SBATCH --partition=open
#SBATCH -o logs/2e_OriginalJordanScriptRecoded.log.out-%a
#SBATCH -e logs/2e_OriginalJordanScriptRecoded.log.err-%a
#SBATCH --array 1-89

# purpose - iake bedfiles of all quartiles from TF/nuc ratio script -> make bedfiles of all quartiles -> runs DNA shape. v5 - 240922; v6 - more code

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/ocean/projects/see180003p/owlang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
WRK=/storage/group/bfp2/default/owl5022-OliviaLang/2023-Krebs_BenzonaseSeq/X_Bulk_Processing
METADATA=../03_Call_Motifs/TF_JASPAR_ENCODE_config.txt
THREADS=4
###

# Dependencies
# - java
# - opencv
# - perl
# - python

set -exo
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/bx

# Load configs
TF=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
JASPAR=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
ENCODE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`

# Inputs and outputs
GENOME=../data/hg38_files/hg38.fa

#set bedfiles
QUARTILE4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/CTCF_NucOccupancy_settings_pipeline_MA1929_1_240910/MA1929_1_final_1000bp_intersected_164bp_category4_1000bp.bed
MEME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240830_TFBS/MEME_files_240904/MA1929.1.meme

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/03_DNAshape_240922/CTCF_DNAshape_MA1929_1_Quartile4_v7_241006

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SMOOTH3=../bin/smoothing_240813.py
EXTRACT=../bin/extract_row_number_240817.py
MASKED=../bin/masked_region_DNAshape_241006.py
MAX_MIN_SCALE=../bin/max_min_scale_v2_241006.py
FORMAT=../bin/format_240922.py
FINAL=../bin/final_240922.py

#set output file names
DNAshape=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HelT=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_HelT.out"}')
MGW=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_MGW.out"}')
PropT=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_PropT.out"}')
Roll=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_AVG_Roll.out"}')
HelT_3=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_HelT_smoothed_3bp.tab"}')
MGW_3=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_MGW_smoothed_3bp.tab"}')
PropT_3=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_PropT_smoothed_3bp.tab"}')
Roll_3=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_Roll_smoothed_3bp.tab"}')
HelT_scale=$(echo "$HelT_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
MGW_scale=$(echo "$MGW_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
PropT_scale=$(echo "$PropT_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
Roll_scale=$(echo "$Roll_3" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale.tab"}')
HelT_final=$(echo "01_""$HelT" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
MGW_final=$(echo "02_""$MGW" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
PropT_final=$(echo "03_""$PropT" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
Roll_final=$(echo "04_""$Roll" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
FILE_final=$(echo "$QUARTILE4" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final_file.tab"}')
NT_count=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NT_count.tab"}')
MASKED_region=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_masked.tab"}')

sampleID=CTCF_DNAshape_Quartile4_v7_241006.slurm

#determine DNA shape with scriptmanager
java -Djava.awt.headless=true -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed --all --avg-composite -o=$DNAshape $GENOME $QUARTILE4
#apply 3 bp smoothing
python $SMOOTH3 $HelT $HelT_3
python $SMOOTH3 $MGW $MGW_3
python $SMOOTH3 $PropT $PropT_3
python $SMOOTH3 $Roll $Roll_3
#extract number of NTs from MEME file
python $EXTRACT $MEME $NT_count
#determine the 5' and 3' boundaries of the motif masked region relative to the center column of tab files at column 256
python $MASKED $NT_count $MASKED_region
#determine max scale for +/- 2bp arbitary units
python $MAX_MIN_SCALE $HelT_3 $MASKED_region $HelT_scale
python $MAX_MIN_SCALE $MGW_3 $MASKED_region $MGW_scale
python $MAX_MIN_SCALE $PropT_3 $MASKED_region $PropT_scale
python $MAX_MIN_SCALE $Roll_3 $MASKED_region $Roll_scale
#change file name of OUT file and swap signs (if most of values are negative)
python $FORMAT $HelT_scale $HelT $HelT_final
python $FORMAT $MGW_scale $MGW $MGW_final
python $FORMAT $PropT_scale $PropT $PropT_final
python $FORMAT $Roll_scale $Roll $Roll_final
#make file of scale and indication if values were swapped
python $FINAL $HelT_scale $MGW_scale $PropT_scale $Roll_scale $FILE_final
