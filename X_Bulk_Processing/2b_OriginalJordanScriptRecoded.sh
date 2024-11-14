#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24GB
#SBATCH --time=1:00:00
#SBATCH --partition=open
#SBATCH -o logs/2b_OriginalJordanScriptRecoded.log.out-%a
#SBATCH -e logs/2b_OriginalJordanScriptRecoded.log.err-%a
#SBATCH --array 1-89

# purpose - intersect bedfile of all sites for a motif with encode-called peaks (Bed narrowPeak bedfile). v8 - intersection (and not intersection) based on Olivia's code. Rest of code from v6. v10 - take ratio, sort, then take quartiles. v11 - updated to hg38

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
BAMFILE=../data/BAM/BNase-seq_50U-10min_merge_hg38.bam		#BAM1
BLACKLIST=../data/hg38_files/ENCFF356LFX_hg38_exclude.bed.gz

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/01_final_TF_pipeline_jobs_240913/CTCF_NucOccupancy_settings_pipeline_MA1929_1_v13_241011

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
JOB=../bin/sum_Col_CDT.pl
JOB_ROW=../bin/sum_Row_CDT.pl
DEDUP=../bin/dedup_coord_by_ID.py
EXTRACT=../bin/extract_row_number_240817.py
MASKED=../bin/masked_region_240817.py
SMOOTH3=../bin/smoothing_240813.py
SMOOTH20=../bin/smoothing_20_240820.py
MAX=../bin/max_position_v3_240818.py
SCALE=../bin/scaling_240814.py
TRANSLATIONAL_sense=../bin/translational_range_sense_v2_240820.py
TRANSLATIONAL_anti=../bin/translational_range_anti_v2_240820.py
TRANSLATIONAL_average=../bin/translational_range_average_240820.py
AUTO=../bin/autocorrelation_of_CDT_v2_240818.py
PERIODICITY=../bin/periodicity_240818.py
ROTATIONAL_sense=../bin/rotational_ratio_sense_v7f_240906.py
MODE_sense=../bin/rotational_sense_mode_v2_240826.py
MODE_sense_substitute=../bin/MODE_sense_substitute_241011.py
PEAKS_shift=../bin/rotational_peaks_shift_v2_241011.py
PEAKS_fill=../bin/rotational_peaks_shifted_columns_v3_240825.py
ROTATIONAL_anti=../bin/rotational_ratio_anti_v7f_240906.py
MODE_anti=../bin/rotational_anti_mode_v2_240826.py
FILTER=../bin/rotational_peaks_filter_240826.py
ROTATIONAL_magnitude=../bin/rotational_magnitude_v2_240826.py
SENSE_count=../bin/rotational_sense_count_240826.py
ANTI_count=../bin/rotational_anti_count_240826.py
CONCAT=../bin/concatenate_v3_241011.py

# ===============================================================================================================================


#set output file names
fileID=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NucOccupancy_settings_pipeline_v13_241011"}')
BEDFILE_a=$(echo "$BEDFILE" | rev | cut -d"/" -f1 | rev | sed 's/_final_1000bp.bed//g' | awk '{print $1}')
ENCODE_BEDFILE_unzipped=$(echo $ENCODE_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1".bed"}')
ENCODE_BEDFILE_shuffled=$(echo $ENCODE_BEDFILE_unzipped | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
BEDFILE_shuffled=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
ENCODE_BEDFILE_1000bp=$(echo $ENCODE_BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BEDFILE_20bp=$(echo $BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_20bp.bed"}')
TARGET_INTERSECT_wDUP=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_intersected_wDUP.bed"}')
TARGET_noINTERSECT_wDUP=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NOTintersected_wDUP.bed"}')
TARGET_Bound=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_intersected.bed"}')
TARGET_noINTERSECT=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NOTintersected.bed"}')
TARGET_Bound_164bp=$(echo $TARGET_Bound  | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_164bp.bed"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1}')
CDT1b=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT1c=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined.cdt"}')
CDT1_sum=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined_sum.tsv"}')
CDT1_sum_noHeader=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined_sum_noHeader.tsv"}')
TSV_ratio=$(echo $TARGET_Bound | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_ratio.TSV"}')
NUMBER=$(echo "$TARGET_Bound" | awk -F. '{print $1"_rowsNumber.tab"}')
NUMBER2=$(echo "$TARGET_Bound" | awk -F. '{print $1"_rowsNumber2.tab"}')
NUMBER3=$(echo "$TARGET_Bound" | awk -F. '{print $1"_rowsNumber3.tab"}')
BEDFILE_category1=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category1.bed"}')
BEDFILE_category2=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category2.bed"}')
BEDFILE_category3=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category3.bed"}')
BEDFILE_category4=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category4.bed"}')
NUMBER_category1=$(echo "$TARGET_Bound_164bp" | awk -F. '{print $1"_category1_rowsNumber.tab"}')
NUMBER_category2=$(echo "$TARGET_Bound_164bp" | awk -F. '{print $1"_category2_rowsNumber.tab"}')
NUMBER_category3=$(echo "$TARGET_Bound_164bp" | awk -F. '{print $1"_category3_rowsNumber.tab"}')
NUMBER_category4=$(echo "$TARGET_Bound_164bp" | awk -F. '{print $1"_category4_rowsNumber.tab"}')
TAB=$(echo "$TARGET_Bound_164bp" | awk -F. '{print $1"_all_rows_values.tab"}')
BEDFILE_category1_1000bp=$(echo $BEDFILE_category1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BEDFILE_category2_1000bp=$(echo $BEDFILE_category2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BEDFILE_category3_1000bp=$(echo $BEDFILE_category3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BEDFILE_category4_1000bp=$(echo $BEDFILE_category4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT2=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT2=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads"}')
CDT2_sense_gz=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT2_anti_gz=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT2_sense=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT2_anti=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
OUT2_sense=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_allReads_sense.tab"}')
OUT2_anti=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_allReads_anti.tab"}')
OUT2_final=$(echo "01_""$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_final.tab"}')
OUT3=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT3=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads"}')
CDT3_sense_gz=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT3_anti_gz=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT3_sense=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT3_anti=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
OUT3_sense=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_allReads_sense.tab"}')
OUT3_anti=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_allReads_anti.tab"}')
OUT3_final=$(echo "02_""$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_final.tab"}')
OUT4=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT4=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads"}')
CDT4_sense_gz=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT4_anti_gz=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT4_sense=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT4_anti=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
OUT4_sense=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_allReads_sense.tab"}')
OUT4_anti=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_allReads_anti.tab"}')
OUT4_final=$(echo "03_""$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_final.tab"}')
OUT5=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT5=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads"}')
CDT5_sense_gz=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT5_anti_gz=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT5_sense=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT5_anti=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
OUT5_sense=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_allReads_sense.tab"}')
OUT5_anti=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_allReads_anti.tab"}')
OUT5_final=$(echo "04_""$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_final.tab"}')
NT_count=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NT_count.tab"}')
MASKED_region=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_masked.tab"}')
category1_sense_smoothed_3=$(echo $OUT2_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category1_anti_smoothed_3=$(echo $OUT2_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category2_sense_smoothed_3=$(echo $OUT3_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category2_anti_smoothed_3=$(echo $OUT3_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category3_sense_smoothed_3=$(echo $OUT4_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category3_anti_smoothed_3=$(echo $OUT4_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category4_sense_smoothed_3=$(echo $OUT5_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category4_anti_smoothed_3=$(echo $OUT5_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_3bp.tab"}')
category1_sense_smoothed_20=$(echo $OUT2_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category1_anti_smoothed_20=$(echo $OUT2_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category2_sense_smoothed_20=$(echo $OUT3_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category2_anti_smoothed_20=$(echo $OUT3_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category3_sense_smoothed_20=$(echo $OUT4_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category3_anti_smoothed_20=$(echo $OUT4_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category4_sense_smoothed_20=$(echo $OUT5_sense | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category4_anti_smoothed_20=$(echo $OUT5_anti | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_smoothed_20bp.tab"}')
category1_sense_max=$(echo $category1_sense_smoothed_20 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_max.tab"}')
category2_sense_max=$(echo $category2_sense_smoothed_20 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_max.tab"}')
category3_sense_max=$(echo $category3_sense_smoothed_20 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_max.tab"}')
category4_sense_max=$(echo $category4_sense_smoothed_20 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_max.tab"}')
all_max_values=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_all_max_values.tab"}')
scale_values=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_scale_values.tab"}')
translational_category1_sense=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category1_sense_translational_setting.tab"}')
translational_category2_sense=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category2_sense_translational_setting.tab"}')
translational_category3_sense=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category3_sense_translational_setting.tab"}')
translational_category4_sense=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category4_sense_translational_setting.tab"}')
translational_category1_anti=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category1_anti_translational_setting.tab"}')
translational_category2_anti=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category2_anti_translational_setting.tab"}')
translational_category3_anti=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category3_anti_translational_setting.tab"}')
translational_category4_anti=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category4_anti_translational_setting.tab"}')
translational_category1=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category1_translational_setting.tab"}')
translational_category2=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category2_translational_setting.tab"}')
translational_category3=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category3_translational_setting.tab"}')
translational_category4=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category4_translational_setting.tab"}')
translational_values=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_translational_values.tab"}')
PERIOD=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_periodicity.tab"}')
correlation_results=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_correlation_results"}')
correlation_results2=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_correlation_results.tsv"}')
rotational_category1_sense=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_sense_rotational_setting.tab"}')
rotational_category1_anti=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_anti_rotational_setting.tab"}')
rotational_category1_sense2=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_sense_rotational_setting_final.tab"}')
rotational_category1_anti2=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_anti_rotational_setting_final.tab"}')
rotational_category1_magnitude=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_rotational_setting_magnitude.tab"}')
rotational_values=$(echo $MEME | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_rotational_values.tab"}')
FINAL=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_final.tab"}')
category1_sense_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category1_sense_smoothed_3_final.tab"}')
category2_sense_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category2_sense_smoothed_3_final.tab"}')
category3_sense_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category3_sense_smoothed_3_final.tab"}')
category4_sense_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category4_sense_smoothed_3_final.tab"}')
category1_anti_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category1_anti_smoothed_3_final.tab"}')
category2_anti_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category2_anti_smoothed_3_final.tab"}')
category3_anti_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category3_anti_smoothed_3_final.tab"}')
category4_anti_smoothed_3_final=$(echo $MEME| rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category4_anti_smoothed_3_final.tab"}')

sampleID=$fileID\.slurm
rm -f $sampleID
$JOBSTATS
#set output
cd $OUTPUT
#unzip files
gunzip -c $ENCODE_BEDFILE > $ENCODE_BEDFILE_unzipped
#shuffle bedfiles
shuf $ENCODE_BEDFILE_unzipped > $ENCODE_BEDFILE_shuffled
shuf $BEDFILE > $BEDFILE_shuffled
#expand befiles
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $ENCODE_BEDFILE_shuffled -o=$OUTPUT/$ENCODE_BEDFILE_1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=20 $BEDFILE_shuffled -o=$OUTPUT/$BEDFILE_20bp
# Intersect peaks with motifs - filter to keep overlap - move ENCODE ChIP value ("signal value") to score col - sort by ID, then score
bedtools intersect -loj -a $BEDFILE_20bp -b $ENCODE_BEDFILE_1000bp | awk '{OFS="\t"}{FS="\t"}{if($8>0) print $1,$2,$3,$4,$13,$6}' | sort -rnk4,5 > $TARGET_INTERSECT_wDUP
bedtools intersect -loj -a $BEDFILE_20bp -b $ENCODE_BEDFILE_1000bp | awk '{OFS="\t"}{FS="\t"}{if($8==-1) print $1,$2,$3,$4,$13,$6}' > $TARGET_noINTERSECT_wDUP
#Deduplicate bound motifs by keeping first instance (larger ENCODE score based on previous command sort)
python $DEDUP -i $TARGET_INTERSECT_wDUP -o $TARGET_Bound
#Deduplicate of unbound motifs does  NOT work as each sites seems to have its own unique 4th column
#get number of rows from intersected bedfile
#expand intersected bedfile
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=164 $TARGET_Bound -o=$OUTPUT/$TARGET_Bound_164bp
#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $TARGET_Bound_164bp $BAM1
#unzip cdt files
gunzip -c $CDT1b > $CDT1c
#no need to scale here as everything is relative to Benzonase data
#sum the number of tags by each row
java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT1_sum -r=1 $CDT1c
#remove header from CDT1_sum file
cat $CDT1_sum | sed '1d' > $CDT1_sum_noHeader
#paste bedfile and CTD1_sum_noHeader, make sure all rows match first, avoid any rows with 0 in TF signal (column 5) or nucleosome occupancy (column 8) then divide encode TF signal to nucleosome occupancy (ratio) in column 7 and sort
paste $TARGET_Bound_164bp $CDT1_sum_noHeader | awk '{OFS="\t"}{FS="\t"}{if ($12=$7 && $5!=0 && $8!=0) print $1,$2,$3,$4,$5,$6,($5/$8)}' | sort -k7,7n > $TSV_ratio
#get number of rows from intersected bedfile
cat $TARGET_Bound | wc -l | awk '{printf "%.f\n", $1 * 0.25}' > $NUMBER
cat $TARGET_Bound | wc -l | awk '{printf "%.f\n", $1 * 0.5}' > $NUMBER2
cat $TARGET_Bound | wc -l | awk '{printf "%.f\n", $1 * 0.75}' > $NUMBER3
#take sorted sites and split into quartiles
#take above motif dedup bedfile that is sorted by TSV.
cat $TSV_ratio | head -$(cat $NUMBER)  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $BEDFILE_category1
cat $TSV_ratio | head -$(cat $NUMBER2)  | tail -$(cat $NUMBER) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $BEDFILE_category2
cat $TSV_ratio | head -$(cat $NUMBER3)  | tail -$(cat $NUMBER) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $BEDFILE_category3
cat $TSV_ratio | tail -$(cat $NUMBER) | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $BEDFILE_category4
#get the number of sites / category
cat $BEDFILE_category1 | wc -l > $NUMBER_category1
cat $BEDFILE_category2 | wc -l > $NUMBER_category2
cat $BEDFILE_category3 | wc -l > $NUMBER_category3
cat $BEDFILE_category4 | wc -l > $NUMBER_category4
#make rows of CSV values of rows / category
echo | awk -v V1="$BEDFILE_a" -v V2=$(cat $BEDFILE_category1 | wc -l) -v V3=$(cat $BEDFILE_category2 | wc -l) -v V4=$(cat $BEDFILE_category3 | wc -l) -v V5=$(cat $BEDFILE_category4 | wc -l) 'BEGIN {print "motif_number\tcategory_1\tcategory_2\tcategory_3\tcategory_4\n"V1"\t"V2"\t"V3"\t"V4"\t"V5}' > $TAB
#expand bedfiles
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category1 -o=$BEDFILE_category1_1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category2 -o=$BEDFILE_category2_1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category3 -o=$BEDFILE_category3_1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category4 -o=$BEDFILE_category4_1000bp
#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT2 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $BEDFILE_category1_1000bp $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT3 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $BEDFILE_category2_1000bp $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT4 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 $BEDFILE_category3_1000bp $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT5 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 $BEDFILE_category4_1000bp $BAM1
#unzip cdt files
gunzip -c $CDT2_sense_gz > $CDT2_sense
gunzip -c $CDT2_anti_gz > $CDT2_anti
gunzip -c $CDT3_sense_gz > $CDT3_sense
gunzip -c $CDT3_anti_gz > $CDT3_anti
gunzip -c $CDT4_sense_gz > $CDT4_sense
gunzip -c $CDT4_anti_gz > $CDT4_anti
gunzip -c $CDT5_sense_gz > $CDT5_sense
gunzip -c $CDT5_anti_gz > $CDT5_anti
#make scaled OUT file for each strand
perl $JOB $CDT2_sense $OUT2_sense
perl $JOB $CDT2_anti $OUT2_anti
perl $JOB $CDT3_sense $OUT3_sense
perl $JOB $CDT3_anti $OUT3_anti
perl $JOB $CDT4_sense $OUT4_sense
perl $JOB $CDT4_anti $OUT4_anti
perl $JOB $CDT5_sense $OUT5_sense
perl $JOB $CDT5_anti $OUT5_anti
#concatenate OUT fles and take lines 1,2,4 to final composite files for each library.
cat $OUT2_sense $OUT2_anti | awk 'NR==1;NR==2;NR==4' > $OUT2_final
cat $OUT3_sense $OUT3_anti | awk 'NR==1;NR==2;NR==4' > $OUT3_final
cat $OUT4_sense $OUT4_anti | awk 'NR==1;NR==2;NR==4' > $OUT4_final
cat $OUT5_sense $OUT5_anti | awk 'NR==1;NR==2;NR==4' > $OUT5_final
#extract number of NTs from MEME file
python $EXTRACT $MEME $NT_count
#determine the 5' and 3' boundaries of the motif masked region relative to the center column of tab files at column 501
python $MASKED $NT_count $MASKED_region
#apply 3 bp smoothing
python $SMOOTH3 $OUT2_sense $category1_sense_smoothed_3
python $SMOOTH3 $OUT2_anti $category1_anti_smoothed_3
python $SMOOTH3 $OUT3_sense $category2_sense_smoothed_3
python $SMOOTH3 $OUT3_anti $category2_anti_smoothed_3
python $SMOOTH3 $OUT4_sense $category3_sense_smoothed_3
python $SMOOTH3 $OUT4_anti $category3_anti_smoothed_3
python $SMOOTH3 $OUT5_sense $category4_sense_smoothed_3
python $SMOOTH3 $OUT5_anti $category4_anti_smoothed_3
#apply 20 bp smoothing
python $SMOOTH20 $OUT2_sense $category1_sense_smoothed_20
python $SMOOTH20 $OUT2_anti $category1_anti_smoothed_20
python $SMOOTH20 $OUT3_sense $category2_sense_smoothed_20
python $SMOOTH20 $OUT3_anti $category2_anti_smoothed_20
python $SMOOTH20 $OUT4_sense $category3_sense_smoothed_20
python $SMOOTH20 $OUT4_anti $category3_anti_smoothed_20
python $SMOOTH20 $OUT5_sense $category4_sense_smoothed_20
python $SMOOTH20 $OUT5_anti $category4_anti_smoothed_20
#get max positions (for later scaling) of sense strand from column 276 (bp -225) - 326 (bp-175) AND determine the bp of the max position. OUTPUT file is name, max value, position of max value
python $MAX $category1_sense_smoothed_20 $category1_sense_max
python $MAX $category2_sense_smoothed_20 $category2_sense_max
python $MAX $category3_sense_smoothed_20 $category3_sense_max
python $MAX $category4_sense_smoothed_20 $category4_sense_max
#combine all above tab files (and remove headers of last 3)
cat $category1_sense_max $category2_sense_max $category3_sense_max $category4_sense_max | awk 'NR==1;NR==2;NR==4;NR==6;NR==8' > $all_max_values
#get scaling value for all categories
python $SCALE $all_max_values $scale_values
#get max range (max-min) from -350 to -150 bp for motif strand
python $TRANSLATIONAL_sense $category1_sense_smoothed_20 $translational_category1_sense
python $TRANSLATIONAL_sense $category2_sense_smoothed_20 $translational_category2_sense
python $TRANSLATIONAL_sense $category3_sense_smoothed_20 $translational_category3_sense
python $TRANSLATIONAL_sense $category4_sense_smoothed_20 $translational_category4_sense
#get max range (max-min) from +150 to +350 bp for opposite strand
python $TRANSLATIONAL_anti $category1_anti_smoothed_20 $translational_category1_anti
python $TRANSLATIONAL_anti $category2_anti_smoothed_20 $translational_category2_anti
python $TRANSLATIONAL_anti $category3_anti_smoothed_20 $translational_category3_anti
python $TRANSLATIONAL_anti $category4_anti_smoothed_20 $translational_category4_anti
#get average of range of translational magnitude from both strands
python $TRANSLATIONAL_average $translational_category1_sense $translational_category1_anti $translational_category1
python $TRANSLATIONAL_average $translational_category2_sense $translational_category2_anti $translational_category2
python $TRANSLATIONAL_average $translational_category3_sense $translational_category3_anti $translational_category3
python $TRANSLATIONAL_average $translational_category4_sense $translational_category4_anti $translational_category4
#combine all above tab files (and add first column of quartile info and header). Output is average of peaks from either strand
cat $translational_category1 $translational_category2 $translational_category3 $translational_category4 | awk 'BEGIN{print "Average_Translational_Magnitude"}1' | awk 'BEGIN{quartile[1]="Quartile"; for(i=2;i<=5;i++) quartile[i]=i-1} {print quartile[NR]"\t"$0} NR>5' > $translational_values
#perform autocorrelation to determine most likely periodicity
python $AUTO -i $category1_sense_smoothed_3 -o $correlation_results
python $PERIODICITY $correlation_results2 $PERIOD
#get significant peaks from category1 sense strand and use those respective bins to call peaks from sense strands of categories 2, 3, and 4
python $ROTATIONAL_sense $category1_sense_smoothed_3 q1_nucleosome_region_sense.tab
python $ROTATIONAL_sense $category2_sense_smoothed_3 q2_nucleosome_region_sense.tab
python $ROTATIONAL_sense $category3_sense_smoothed_3 q3_nucleosome_region_sense.tab
python $ROTATIONAL_sense $category4_sense_smoothed_3 q4_nucleosome_region_sense.tab
#get concatenated list of all unique, significant peaks, THEN get their the mode of their max position (bp)
cat q1_nucleosome_region_sense.tab q2_nucleosome_region_sense.tab q3_nucleosome_region_sense.tab q4_nucleosome_region_sense.tab > significant_peaks_sense.tab
python $MODE_sense significant_peaks_sense.tab $MASKED_region significant_peaks_sense_mode.tab
#determine unique set of 'significant' peaks from motif strand and do unique sort
cat q1_nucleosome_region_sense.tab q2_nucleosome_region_sense.tab q3_nucleosome_region_sense.tab q4_nucleosome_region_sense.tab | cut -f1,2 | sort -k1,1 | uniq > output_filtered_nucleosome_sense.tab
#Sense mode needs substituted to match 'opposite strand phase (0-9) then shift by doing by 5' or 3' by 'mode-5' with +=shift 5', -=shift3'
python $MODE_sense_substitute significant_peaks_sense_mode.tab significant_peaks_sense_mode_substituted.tab
#sort unique, significant peaks by the above substituted mode'
python $PEAKS_shift output_filtered_nucleosome_sense.tab significant_peaks_sense_mode_substituted.tab shifted_columns_sense.tab
#take shifted, significant bins and fill out range for each bin
python $PEAKS_fill $category1_sense_smoothed_3 shifted_columns_sense.tab category1_sense_smoothed_3_full.tab
python $PEAKS_fill $category2_sense_smoothed_3 shifted_columns_sense.tab category2_sense_smoothed_3_full.tab
python $PEAKS_fill $category3_sense_smoothed_3 shifted_columns_sense.tab category3_sense_smoothed_3_full.tab
python $PEAKS_fill $category4_sense_smoothed_3 shifted_columns_sense.tab category4_sense_smoothed_3_full.tab
#repeat above for opposite strand
#get significant peaks from category1 anti strand and use those respective bins to call peaks from sense strands of categories 2, 3, and 4
python $ROTATIONAL_anti $category1_anti_smoothed_3 q1_nucleosome_region_anti.tab
python $ROTATIONAL_anti $category2_anti_smoothed_3 q2_nucleosome_region_anti.tab
python $ROTATIONAL_anti $category3_anti_smoothed_3 q3_nucleosome_region_anti.tab
python $ROTATIONAL_anti $category4_anti_smoothed_3 q4_nucleosome_region_anti.tab
#get concatenated list of all unique, significant peaks, THEN get their the mode of their max position (bp)
cat q1_nucleosome_region_anti.tab q2_nucleosome_region_anti.tab q3_nucleosome_region_anti.tab q4_nucleosome_region_anti.tab > significant_peaks_anti.tab
python $MODE_anti significant_peaks_anti.tab $MASKED_region significant_peaks_anti_mode.tab
#determine unique set of 'significant' peaks from opposite strand and do unique sort
cat q1_nucleosome_region_anti.tab q2_nucleosome_region_anti.tab q3_nucleosome_region_anti.tab q4_nucleosome_region_anti.tab | cut -f1,2 | sort -k1,1 | uniq > output_filtered_nucleosome_anti.tab
#sort unique, significant peaks by the above mode; shift is by adding mode/2 to shift 3' to the motif
python $PEAKS_shift output_filtered_nucleosome_anti.tab significant_peaks_anti_mode.tab shifted_columns_anti.tab
#take shifted, significant bins and fill out range for each bin
python $PEAKS_fill $category1_anti_smoothed_3 shifted_columns_anti.tab category1_anti_smoothed_3_full.tab
python $PEAKS_fill $category2_anti_smoothed_3 shifted_columns_anti.tab category2_anti_smoothed_3_full.tab
python $PEAKS_fill $category3_anti_smoothed_3 shifted_columns_anti.tab category3_anti_smoothed_3_full.tab
python $PEAKS_fill $category4_anti_smoothed_3 shifted_columns_anti.tab category4_anti_smoothed_3_full.tab
#remove rows whose max value (column 8) is within the masked motif region
python $FILTER category1_sense_smoothed_3_full.tab $MASKED_region $category1_sense_smoothed_3_final
python $FILTER category2_sense_smoothed_3_full.tab $MASKED_region $category2_sense_smoothed_3_final
python $FILTER category3_sense_smoothed_3_full.tab $MASKED_region $category3_sense_smoothed_3_final
python $FILTER category4_sense_smoothed_3_full.tab $MASKED_region $category4_sense_smoothed_3_final
python $FILTER category1_anti_smoothed_3_full.tab $MASKED_region $category1_anti_smoothed_3_final
python $FILTER category2_anti_smoothed_3_full.tab $MASKED_region $category2_anti_smoothed_3_final
python $FILTER category3_anti_smoothed_3_full.tab $MASKED_region $category3_anti_smoothed_3_final
python $FILTER category4_anti_smoothed_3_full.tab $MASKED_region $category4_anti_smoothed_3_final
#get average of all peaks 5' to motif (motif strand) and 5' to motif (opposite strand)
#calculate average range (magnitude of rotational setting / category) of all peaks on either strand
python $ROTATIONAL_magnitude $category1_sense_smoothed_3_final $category1_anti_smoothed_3_final category1_rotational_magnitude.tab
python $ROTATIONAL_magnitude $category2_sense_smoothed_3_final $category2_anti_smoothed_3_final category2_rotational_magnitude.tab
python $ROTATIONAL_magnitude $category3_sense_smoothed_3_final $category3_anti_smoothed_3_final category3_rotational_magnitude.tab
python $ROTATIONAL_magnitude $category4_sense_smoothed_3_final $category4_anti_smoothed_3_final category4_rotational_magnitude.tab
#combine all above tab files (and add first column of quartile info and header). Output is average of peaks from either strand
cat category1_rotational_magnitude.tab category2_rotational_magnitude.tab category3_rotational_magnitude.tab category4_rotational_magnitude.tab | awk 'NR % 2 == 0' | awk 'BEGIN{print "Average_Rotational_Magnitude"}1' | awk 'BEGIN{quartile[1]="Quartile"; for(i=2;i<=5;i++) quartile[i]=i-1} {print quartile[NR]"\t"$0} NR>5' > $rotational_values
#calculate the # of unique, significant rotational peaks 5' to motif (from motif set of peaks).
python $SENSE_count $category1_sense_smoothed_3_final SENSE_count.tab
#calculate the # of unique, significant rotational peaks 3' to motif (from motif set of peaks).
python $ANTI_count $category1_anti_smoothed_3_final ANTI_count.tab
#make final file with all key information for this TF
python $CONCAT $scale_values $PERIOD $translational_values significant_peaks_sense_mode_substituted.tab significant_peaks_anti_mode.tab $rotational_values SENSE_count.tab ANTI_count.tab $FINAL
