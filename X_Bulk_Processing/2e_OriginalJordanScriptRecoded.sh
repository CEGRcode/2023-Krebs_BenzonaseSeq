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
TARGET=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $1}'`
JASPAR=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $2}'`
ENCODE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $METADATA | awk '{print $3}'`

# Inputs and outputs
MOTIF=../data/RefPT-JASPAR
GENOME=../data/hg38_files/hg38.fa

RUNID=${TARGET}_${JASPAR}
BEDFILE=$MOTIF/1000bp/${RUNID}_SORT-TFnucRatio_GROUP-Quartile4_1000bp.bed
BED=`basename $BEDFILE ".bed"`

RUNID=${TARGET}_${JASPAR}
ODIR=Library/DNAShape_${RUNID}
[ -d $ODIR ] || mkdir -p $ODIR
[ -d logs ] || mkdir logs

BASE=$ODIR/$BED
MASKED_region=${BASE}_masked.tab

#set bedfiles
QUARTILE4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/CTCF_NucOccupancy_settings_pipeline_MA1929_1_240910/MA1929_1_final_1000bp_intersected_164bp_category4_1000bp.bed
MEME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240830_TFBS/MEME_files_240904/MA1929.1.meme

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/03_DNAshape_240922/CTCF_DNAshape_MA1929_1_Quartile4_v7_241006

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
SMOOTH=../bin/smoothing_parameterize.py
EXTRACT=../bin/extract_row_number_240817.py
MASKED=../bin/masked_region_DNAshape_241006.py
MAX_MIN_SCALE=../bin/max_min_scale_v2_241006.py
FORMAT=../bin/format_240922.py
FINAL=../bin/final_240922.py

#extract number of NTs from MEME file
python $EXTRACT $MEME ${BASE}_NT_count.tab

#determine the 5' and 3' boundaries of the motif masked region relative to the center column of tab files at column 256
python $MASKED ${BASE}_NT_count.tab $MASKED_region

#determine DNA shape with scriptmanager
java -Djava.awt.headless=true -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed --all --composite -o $BASE $GENOME $BEDFILE

#apply 3 bp smoothing
python $SMOOTH 3 ${BASE}_AVG_HelT.out ${BASE}_HelT_smooth3.tab
python $SMOOTH 3 ${BASE}_AVG_MGW.out ${BASE}_MGW_smooth3.tab
python $SMOOTH 3 ${BASE}_AVG_PropT.out ${BASE}_PropT_smooth3.tab
python $SMOOTH 3 ${BASE}_AVG_Roll.out ${BASE}_Roll_smooth3.tab

#determine max scale for +/- 2bp arbitary units
python $MAX_MIN_SCALE ${BASE}_HelT_smooth3.tab $MASKED_region ${BASE}_HelT_smooth3_scale.tab
python $MAX_MIN_SCALE ${BASE}_MGW_smooth3.tab $MASKED_region ${BASE}_MGW_smooth3_scale.tab
python $MAX_MIN_SCALE ${BASE}_PropT_smooth3.tab $MASKED_region ${BASE}_PropT_smooth3_scale.tab
python $MAX_MIN_SCALE ${BASE}_Roll_smooth3.tab $MASKED_region ${BASE}_Roll_smooth3_scale.tab

#change file name of OUT file and swap signs (if most of values are negative)
python $FORMAT ${BASE}_HelT_smooth3_scale.tab ${BASE}_AVG_HelT.out $ODIR/01_${BED}_AVG_HelT_final.tab
python $FORMAT ${BASE}_MGW_smooth3_scale.tab ${BASE}_AVG_MGW.out $ODIR/02_${BED}_AVG_MGW_final.tab
python $FORMAT ${BASE}_PropT_smooth3_scale.tab ${BASE}_AVG_PropT.out $ODIR/03_${BED}_AVG_PropT_final.tab
python $FORMAT ${BASE}_Roll_smooth3_scale.tab ${BASE}_AVG_Roll.out $ODIR/04_${BED}_AVG_Roll_final.tab

#make file of scale and indication if values were swapped
python $FINAL ${BASE}_HelT_smooth3_scale.tab \
			${BASE}_MGW_smooth3_scale.tab \
			${BASE}_PropT_smooth3_scale.tab \
			${BASE}_Roll_smooth3_scale.tab \
			${BASE}_final_file.tab
