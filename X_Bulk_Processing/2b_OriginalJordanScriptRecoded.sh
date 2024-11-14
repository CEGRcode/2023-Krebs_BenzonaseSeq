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
echo "$JOBSTATS" >> $sampleID
echo "#set output" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#unzip files" >> $sampleID
echo "gunzip -c $ENCODE_BEDFILE > $ENCODE_BEDFILE_unzipped" >> $sampleID
echo "#shuffle bedfiles" >> $sampleID
echo "shuf $ENCODE_BEDFILE_unzipped > $ENCODE_BEDFILE_shuffled" >> $sampleID
echo "shuf $BEDFILE > $BEDFILE_shuffled" >> $sampleID
echo "#expand befiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $ENCODE_BEDFILE_shuffled -o=$OUTPUT/$ENCODE_BEDFILE_1000bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=20 $BEDFILE_shuffled -o=$OUTPUT/$BEDFILE_20bp" >> $sampleID
echo "# Intersect peaks with motifs - filter to keep overlap - move ENCODE ChIP value ("signal value") to score col - sort by ID, then score" >> $sampleID
echo "bedtools intersect -loj -a $BEDFILE_20bp -b $ENCODE_BEDFILE_1000bp | awk '{OFS=\"\t\"}{FS=\"\t\"}{if(\$8>0) print \$1,\$2,\$3,\$4,\$13,\$6}' | sort -rnk4,5 > $TARGET_INTERSECT_wDUP" >> $sampleID
echo "bedtools intersect -loj -a $BEDFILE_20bp -b $ENCODE_BEDFILE_1000bp | awk '{OFS=\"\t\"}{FS=\"\t\"}{if(\$8==-1) print \$1,\$2,\$3,\$4,\$13,\$6}' > $TARGET_noINTERSECT_wDUP" >> $sampleID
echo "#Deduplicate bound motifs by keeping first instance (larger ENCODE score based on previous command sort)" >> $sampleID
echo "python $DEDUP -i $TARGET_INTERSECT_wDUP -o $TARGET_Bound" >> $sampleID
echo "#Deduplicate of unbound motifs does  NOT work as each sites seems to have its own unique 4th column" >> $sampleID
echo "#get number of rows from intersected bedfile" >> $sampleID
echo "#expand intersected bedfile" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=164 $TARGET_Bound -o=$OUTPUT/$TARGET_Bound_164bp" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $TARGET_Bound_164bp $BAM1" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1b > $CDT1c" >> $sampleID
echo "#no need to scale here as everything is relative to Benzonase data" >> $sampleID
echo "#sum the number of tags by each row" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT1_sum -r=1 $CDT1c" >> $sampleID
echo "#remove header from CDT1_sum file" >> $sampleID
echo "cat $CDT1_sum | sed '1d' > $CDT1_sum_noHeader" >> $sampleID
echo "#paste bedfile and CTD1_sum_noHeader, make sure all rows match first, avoid any rows with 0 in TF signal (column 5) or nucleosome occupancy (column 8) then divide encode TF signal to nucleosome occupancy (ratio) in column 7 and sort" >> $sampleID
echo "paste $TARGET_Bound_164bp $CDT1_sum_noHeader | awk '{OFS=\"\t\"}{FS=\"\t\"}{if (\$12=\$7 && \$5!=0 && \$8!=0) print \$1,\$2,\$3,\$4,\$5,\$6,(\$5/\$8)}' | sort -k7,7n > $TSV_ratio" >> $sampleID
echo "#get number of rows from intersected bedfile" >> $sampleID
echo "cat $TARGET_Bound | wc -l | awk '{printf \"%.f\\n\", \$1 * 0.25}' > $NUMBER" >> $sampleID
echo "cat $TARGET_Bound | wc -l | awk '{printf \"%.f\\n\", \$1 * 0.5}' > $NUMBER2" >> $sampleID
echo "cat $TARGET_Bound | wc -l | awk '{printf \"%.f\\n\", \$1 * 0.75}' > $NUMBER3" >> $sampleID
echo "#take sorted sites and split into quartiles" >> $sampleID
echo "#take above motif dedup bedfile that is sorted by TSV." >> $sampleID
echo "cat $TSV_ratio | head -\$(cat $NUMBER)  | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category1" >> $sampleID
echo "cat $TSV_ratio | head -\$(cat $NUMBER2)  | tail -\$(cat $NUMBER) | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category2" >> $sampleID
echo "cat $TSV_ratio | head -\$(cat $NUMBER3)  | tail -\$(cat $NUMBER) | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category3" >> $sampleID
echo "cat $TSV_ratio | tail -\$(cat $NUMBER) | awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category4" >> $sampleID
echo "#get the number of sites / category" >> $sampleID
echo "cat $BEDFILE_category1 | wc -l > $NUMBER_category1" >> $sampleID
echo "cat $BEDFILE_category2 | wc -l > $NUMBER_category2" >> $sampleID
echo "cat $BEDFILE_category3 | wc -l > $NUMBER_category3" >> $sampleID
echo "cat $BEDFILE_category4 | wc -l > $NUMBER_category4" >> $sampleID
echo "#make rows of CSV values of rows / category" >> $sampleID
echo "echo | awk -v V1="$BEDFILE_a" -v V2=\$(cat $BEDFILE_category1 | wc -l) -v V3=\$(cat $BEDFILE_category2 | wc -l) -v V4=\$(cat $BEDFILE_category3 | wc -l) -v V5=\$(cat $BEDFILE_category4 | wc -l) 'BEGIN {print \"motif_number\"\"\t\"\"category_1\"\"\t\"\"category_2\"\"\t\"\"category_3\"\"\t\"\"category_4\"\"\n\"V1\"\t\"V2\"\t\"V3\"\t\"V4\"\t\"V5}' > $TAB" >> $sampleID
echo "#expand bedfiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category1 -o=$BEDFILE_category1_1000bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category2 -o=$BEDFILE_category2_1000bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category3 -o=$BEDFILE_category3_1000bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $BEDFILE_category4 -o=$BEDFILE_category4_1000bp" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT2 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $BEDFILE_category1_1000bp $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT3 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $BEDFILE_category2_1000bp $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT4 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 $BEDFILE_category3_1000bp $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT5 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 $BEDFILE_category4_1000bp $BAM1" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT2_sense_gz > $CDT2_sense" >> $sampleID
echo "gunzip -c $CDT2_anti_gz > $CDT2_anti" >> $sampleID
echo "gunzip -c $CDT3_sense_gz > $CDT3_sense" >> $sampleID
echo "gunzip -c $CDT3_anti_gz > $CDT3_anti" >> $sampleID
echo "gunzip -c $CDT4_sense_gz > $CDT4_sense" >> $sampleID
echo "gunzip -c $CDT4_anti_gz > $CDT4_anti" >> $sampleID
echo "gunzip -c $CDT5_sense_gz > $CDT5_sense" >> $sampleID
echo "gunzip -c $CDT5_anti_gz > $CDT5_anti" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT2_sense $OUT2_sense" >> $sampleID
echo "perl $JOB $CDT2_anti $OUT2_anti" >> $sampleID
echo "perl $JOB $CDT3_sense $OUT3_sense" >> $sampleID
echo "perl $JOB $CDT3_anti $OUT3_anti" >> $sampleID
echo "perl $JOB $CDT4_sense $OUT4_sense" >> $sampleID
echo "perl $JOB $CDT4_anti $OUT4_anti" >> $sampleID
echo "perl $JOB $CDT5_sense $OUT5_sense" >> $sampleID
echo "perl $JOB $CDT5_anti $OUT5_anti" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library." >> $sampleID
echo "cat $OUT2_sense $OUT2_anti | awk 'NR==1;NR==2;NR==4' > $OUT2_final" >> $sampleID
echo "cat $OUT3_sense $OUT3_anti | awk 'NR==1;NR==2;NR==4' > $OUT3_final" >> $sampleID
echo "cat $OUT4_sense $OUT4_anti | awk 'NR==1;NR==2;NR==4' > $OUT4_final" >> $sampleID
echo "cat $OUT5_sense $OUT5_anti | awk 'NR==1;NR==2;NR==4' > $OUT5_final" >> $sampleID
echo "#extract number of NTs from MEME file" >> $sampleID
echo "python3 $EXTRACT $MEME $NT_count" >> $sampleID
echo "#determine the 5' and 3' boundaries of the motif masked region relative to the center column of tab files at column 501" >> $sampleID
echo "python3 $MASKED $NT_count $MASKED_region" >> $sampleID
echo "#apply 3 bp smoothing" >> $sampleID
echo "python3 $SMOOTH3 $OUT2_sense $category1_sense_smoothed_3" >> $sampleID
echo "python3 $SMOOTH3 $OUT2_anti $category1_anti_smoothed_3" >> $sampleID
echo "python3 $SMOOTH3 $OUT3_sense $category2_sense_smoothed_3" >> $sampleID
echo "python3 $SMOOTH3 $OUT3_anti $category2_anti_smoothed_3" >> $sampleID
echo "python3 $SMOOTH3 $OUT4_sense $category3_sense_smoothed_3" >> $sampleID
echo "python3 $SMOOTH3 $OUT4_anti $category3_anti_smoothed_3" >> $sampleID
echo "python3 $SMOOTH3 $OUT5_sense $category4_sense_smoothed_3" >> $sampleID
echo "python3 $SMOOTH3 $OUT5_anti $category4_anti_smoothed_3" >> $sampleID
echo "#apply 20 bp smoothing" >> $sampleID
echo "python3 $SMOOTH20 $OUT2_sense $category1_sense_smoothed_20" >> $sampleID
echo "python3 $SMOOTH20 $OUT2_anti $category1_anti_smoothed_20" >> $sampleID
echo "python3 $SMOOTH20 $OUT3_sense $category2_sense_smoothed_20" >> $sampleID
echo "python3 $SMOOTH20 $OUT3_anti $category2_anti_smoothed_20" >> $sampleID
echo "python3 $SMOOTH20 $OUT4_sense $category3_sense_smoothed_20" >> $sampleID
echo "python3 $SMOOTH20 $OUT4_anti $category3_anti_smoothed_20" >> $sampleID
echo "python3 $SMOOTH20 $OUT5_sense $category4_sense_smoothed_20" >> $sampleID
echo "python3 $SMOOTH20 $OUT5_anti $category4_anti_smoothed_20" >> $sampleID
echo "#get max positions (for later scaling) of sense strand from column 276 (bp -225) - 326 (bp-175) AND determine the bp of the max position. OUTPUT file is name, max value, position of max value" >> $sampleID
echo "python3 $MAX $category1_sense_smoothed_20 $category1_sense_max" >> $sampleID
echo "python3 $MAX $category2_sense_smoothed_20 $category2_sense_max" >> $sampleID
echo "python3 $MAX $category3_sense_smoothed_20 $category3_sense_max" >> $sampleID
echo "python3 $MAX $category4_sense_smoothed_20 $category4_sense_max" >> $sampleID
echo "#combine all above tab files (and remove headers of last 3)" >> $sampleID
echo "cat $category1_sense_max $category2_sense_max $category3_sense_max $category4_sense_max | awk 'NR==1;NR==2;NR==4;NR==6;NR==8' > $all_max_values" >> $sampleID
echo "#get scaling value for all categories" >> $sampleID
echo "python3 $SCALE $all_max_values $scale_values" >> $sampleID
echo "#get max range (max-min) from -350 to -150 bp for motif strand" >> $sampleID
echo "python3 $TRANSLATIONAL_sense $category1_sense_smoothed_20 $translational_category1_sense" >> $sampleID
echo "python3 $TRANSLATIONAL_sense $category2_sense_smoothed_20 $translational_category2_sense" >> $sampleID
echo "python3 $TRANSLATIONAL_sense $category3_sense_smoothed_20 $translational_category3_sense" >> $sampleID
echo "python3 $TRANSLATIONAL_sense $category4_sense_smoothed_20 $translational_category4_sense" >> $sampleID
echo "#get max range (max-min) from +150 to +350 bp for opposite strand" >> $sampleID
echo "python3 $TRANSLATIONAL_anti $category1_anti_smoothed_20 $translational_category1_anti" >> $sampleID
echo "python3 $TRANSLATIONAL_anti $category2_anti_smoothed_20 $translational_category2_anti" >> $sampleID
echo "python3 $TRANSLATIONAL_anti $category3_anti_smoothed_20 $translational_category3_anti" >> $sampleID
echo "python3 $TRANSLATIONAL_anti $category4_anti_smoothed_20 $translational_category4_anti" >> $sampleID
echo "#get average of range of translational magnitude from both strands" >> $sampleID
echo "python3 $TRANSLATIONAL_average $translational_category1_sense $translational_category1_anti $translational_category1" >> $sampleID
echo "python3 $TRANSLATIONAL_average $translational_category2_sense $translational_category2_anti $translational_category2" >> $sampleID
echo "python3 $TRANSLATIONAL_average $translational_category3_sense $translational_category3_anti $translational_category3" >> $sampleID
echo "python3 $TRANSLATIONAL_average $translational_category4_sense $translational_category4_anti $translational_category4" >> $sampleID
echo "#combine all above tab files (and add first column of quartile info and header). Output is average of peaks from either strand" >> $sampleID
echo "cat $translational_category1 $translational_category2 $translational_category3 $translational_category4 | awk 'BEGIN{print \"Average_Translational_Magnitude\"}1' | awk 'BEGIN{quartile[1]=\"Quartile\"; for(i=2;i<=5;i++) quartile[i]=i-1} {print quartile[NR]\"\t\"\$0} NR>5' > $translational_values" >> $sampleID
echo "#perform autocorrelation to determine most likely periodicity" >> $sampleID
echo "python3 $AUTO -i $category1_sense_smoothed_3 -o $correlation_results" >> $sampleID
echo "python $PERIODICITY $correlation_results2 $PERIOD" >> $sampleID
echo "#get significant peaks from category1 sense strand and use those respective bins to call peaks from sense strands of categories 2, 3, and 4" >> $sampleID
echo "python3 $ROTATIONAL_sense $category1_sense_smoothed_3 q1_nucleosome_region_sense.tab" >> $sampleID
echo "python3 $ROTATIONAL_sense $category2_sense_smoothed_3 q2_nucleosome_region_sense.tab" >> $sampleID
echo "python3 $ROTATIONAL_sense $category3_sense_smoothed_3 q3_nucleosome_region_sense.tab" >> $sampleID
echo "python3 $ROTATIONAL_sense $category4_sense_smoothed_3 q4_nucleosome_region_sense.tab" >> $sampleID
echo "#get concatenated list of all unique, significant peaks, THEN get their the mode of their max position (bp)" >> $sampleID
echo "cat q1_nucleosome_region_sense.tab q2_nucleosome_region_sense.tab q3_nucleosome_region_sense.tab q4_nucleosome_region_sense.tab > significant_peaks_sense.tab" >> $sampleID
echo "python3 $MODE_sense significant_peaks_sense.tab $MASKED_region significant_peaks_sense_mode.tab" >> $sampleID
echo "#determine unique set of 'significant' peaks from motif strand and do unique sort" >> $sampleID
echo "cat q1_nucleosome_region_sense.tab q2_nucleosome_region_sense.tab q3_nucleosome_region_sense.tab q4_nucleosome_region_sense.tab | cut -f1,2 | sort -k1,1 | uniq > output_filtered_nucleosome_sense.tab" >> $sampleID
echo "#Sense mode needs substituted to match 'opposite strand phase (0-9) then shift by doing by 5' or 3' by 'mode-5' with +=shift 5', -=shift3'" >> $sampleID
echo "python3 $MODE_sense_substitute significant_peaks_sense_mode.tab significant_peaks_sense_mode_substituted.tab" >> $sampleID
echo "#sort unique, significant peaks by the above substituted mode'" >> $sampleID
echo "python3 $PEAKS_shift output_filtered_nucleosome_sense.tab significant_peaks_sense_mode_substituted.tab shifted_columns_sense.tab" >> $sampleID
echo "#take shifted, significant bins and fill out range for each bin" >> $sampleID
echo "python3 $PEAKS_fill $category1_sense_smoothed_3 shifted_columns_sense.tab category1_sense_smoothed_3_full.tab" >> $sampleID
echo "python3 $PEAKS_fill $category2_sense_smoothed_3 shifted_columns_sense.tab category2_sense_smoothed_3_full.tab" >> $sampleID
echo "python3 $PEAKS_fill $category3_sense_smoothed_3 shifted_columns_sense.tab category3_sense_smoothed_3_full.tab" >> $sampleID
echo "python3 $PEAKS_fill $category4_sense_smoothed_3 shifted_columns_sense.tab category4_sense_smoothed_3_full.tab" >> $sampleID
echo "#repeat above for opposite strand" >> $sampleID
echo "#get significant peaks from category1 anti strand and use those respective bins to call peaks from sense strands of categories 2, 3, and 4" >> $sampleID
echo "python3 $ROTATIONAL_anti $category1_anti_smoothed_3 q1_nucleosome_region_anti.tab" >> $sampleID
echo "python3 $ROTATIONAL_anti $category2_anti_smoothed_3 q2_nucleosome_region_anti.tab" >> $sampleID
echo "python3 $ROTATIONAL_anti $category3_anti_smoothed_3 q3_nucleosome_region_anti.tab" >> $sampleID
echo "python3 $ROTATIONAL_anti $category4_anti_smoothed_3 q4_nucleosome_region_anti.tab" >> $sampleID
echo "#get concatenated list of all unique, significant peaks, THEN get their the mode of their max position (bp)" >> $sampleID
echo "cat q1_nucleosome_region_anti.tab q2_nucleosome_region_anti.tab q3_nucleosome_region_anti.tab q4_nucleosome_region_anti.tab > significant_peaks_anti.tab" >> $sampleID
echo "python3 $MODE_anti significant_peaks_anti.tab $MASKED_region significant_peaks_anti_mode.tab" >> $sampleID
echo "#determine unique set of 'significant' peaks from opposite strand and do unique sort" >> $sampleID
echo "cat q1_nucleosome_region_anti.tab q2_nucleosome_region_anti.tab q3_nucleosome_region_anti.tab q4_nucleosome_region_anti.tab | cut -f1,2 | sort -k1,1 | uniq > output_filtered_nucleosome_anti.tab" >> $sampleID
echo "#sort unique, significant peaks by the above mode; shift is by adding mode/2 to shift 3' to the motif" >> $sampleID
echo "python3 $PEAKS_shift output_filtered_nucleosome_anti.tab significant_peaks_anti_mode.tab shifted_columns_anti.tab" >> $sampleID
echo "#take shifted, significant bins and fill out range for each bin" >> $sampleID
echo "python3 $PEAKS_fill $category1_anti_smoothed_3 shifted_columns_anti.tab category1_anti_smoothed_3_full.tab" >> $sampleID
echo "python3 $PEAKS_fill $category2_anti_smoothed_3 shifted_columns_anti.tab category2_anti_smoothed_3_full.tab" >> $sampleID
echo "python3 $PEAKS_fill $category3_anti_smoothed_3 shifted_columns_anti.tab category3_anti_smoothed_3_full.tab" >> $sampleID
echo "python3 $PEAKS_fill $category4_anti_smoothed_3 shifted_columns_anti.tab category4_anti_smoothed_3_full.tab" >> $sampleID
echo "#remove rows whose max value (column 8) is within the masked motif region" >> $sampleID
echo "python3 $FILTER category1_sense_smoothed_3_full.tab $MASKED_region $category1_sense_smoothed_3_final" >> $sampleID
echo "python3 $FILTER category2_sense_smoothed_3_full.tab $MASKED_region $category2_sense_smoothed_3_final" >> $sampleID
echo "python3 $FILTER category3_sense_smoothed_3_full.tab $MASKED_region $category3_sense_smoothed_3_final" >> $sampleID
echo "python3 $FILTER category4_sense_smoothed_3_full.tab $MASKED_region $category4_sense_smoothed_3_final" >> $sampleID
echo "python3 $FILTER category1_anti_smoothed_3_full.tab $MASKED_region $category1_anti_smoothed_3_final" >> $sampleID
echo "python3 $FILTER category2_anti_smoothed_3_full.tab $MASKED_region $category2_anti_smoothed_3_final" >> $sampleID
echo "python3 $FILTER category3_anti_smoothed_3_full.tab $MASKED_region $category3_anti_smoothed_3_final" >> $sampleID
echo "python3 $FILTER category4_anti_smoothed_3_full.tab $MASKED_region $category4_anti_smoothed_3_final" >> $sampleID
echo "#get average of all peaks 5' to motif (motif strand) and 5' to motif (opposite strand)" >> $sampleID
echo "#calculate average range (magnitude of rotational setting / category) of all peaks on either strand" >> $sampleID
echo "python3 $ROTATIONAL_magnitude $category1_sense_smoothed_3_final $category1_anti_smoothed_3_final category1_rotational_magnitude.tab" >> $sampleID
echo "python3 $ROTATIONAL_magnitude $category2_sense_smoothed_3_final $category2_anti_smoothed_3_final category2_rotational_magnitude.tab" >> $sampleID
echo "python3 $ROTATIONAL_magnitude $category3_sense_smoothed_3_final $category3_anti_smoothed_3_final category3_rotational_magnitude.tab" >> $sampleID
echo "python3 $ROTATIONAL_magnitude $category4_sense_smoothed_3_final $category4_anti_smoothed_3_final category4_rotational_magnitude.tab" >> $sampleID
echo "#combine all above tab files (and add first column of quartile info and header). Output is average of peaks from either strand" >> $sampleID
echo "cat category1_rotational_magnitude.tab category2_rotational_magnitude.tab category3_rotational_magnitude.tab category4_rotational_magnitude.tab | awk 'NR % 2 == 0' | awk 'BEGIN{print \"Average_Rotational_Magnitude\"}1' | awk 'BEGIN{quartile[1]=\"Quartile\"; for(i=2;i<=5;i++) quartile[i]=i-1} {print quartile[NR]\"\t\"\$0} NR>5' > $rotational_values" >> $sampleID
echo "#calculate the # of unique, significant rotational peaks 5' to motif (from motif set of peaks)." >> $sampleID
echo "python $SENSE_count $category1_sense_smoothed_3_final SENSE_count.tab" >> $sampleID
echo "#calculate the # of unique, significant rotational peaks 3' to motif (from motif set of peaks)." >> $sampleID
echo "python $ANTI_count $category1_anti_smoothed_3_final ANTI_count.tab" >> $sampleID
echo "#make final file with all key information for this TF" >> $sampleID
echo "python3 $CONCAT $scale_values $PERIOD $translational_values significant_peaks_sense_mode_substituted.tab significant_peaks_anti_mode.tab $rotational_values SENSE_count.tab ANTI_count.tab $FINAL" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $ENCODE_BEDFILE_unzipped" >> $sampleID
echo "rm $ENCODE_BEDFILE_shuffled" >> $sampleID
echo "rm $BEDFILE_shuffled" >> $sampleID
echo "rm $ENCODE_BEDFILE_1000bp" >> $sampleID
echo "rm $BEDFILE_20bp" >> $sampleID
echo "rm $TARGET_INTERSECT_wDUP" >> $sampleID
echo "rm $TARGET_noINTERSECT_wDUP" >> $sampleID
echo "rm $TARGET_Bound" >> $sampleID
echo "rm $TARGET_noINTERSECT" >> $sampleID
echo "rm $TARGET_Bound_164bp" >> $sampleID
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1" >> $sampleID
echo "rm $CDT1b" >> $sampleID
echo "rm $CDT1c" >> $sampleID
echo "rm $CDT1_sum" >> $sampleID
echo "rm $CDT1_sum_noHeader" >> $sampleID
echo "rm $NUMBER" >> $sampleID
echo "rm $NUMBER2" >> $sampleID
echo "rm $NUMBER3" >> $sampleID
echo "rm $BEDFILE_category1" >> $sampleID
echo "rm $BEDFILE_category2" >> $sampleID
echo "rm $BEDFILE_category3" >> $sampleID
echo "rm $BEDFILE_category4" >> $sampleID
echo "rm $NUMBER_category1" >> $sampleID
echo "rm $NUMBER_category2" >> $sampleID
echo "rm $NUMBER_category3" >> $sampleID
echo "rm $NUMBER_category4" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2" >> $sampleID
echo "rm $CDT2_sense_gz" >> $sampleID
echo "rm $CDT2_anti_gz" >> $sampleID
echo "rm $CDT2_sense" >> $sampleID
echo "rm $CDT2_anti" >> $sampleID
echo "rm $OUT3" >> $sampleID
echo "rm $CDT3" >> $sampleID
echo "rm $CDT3_sense_gz" >> $sampleID
echo "rm $CDT3_anti_gz" >> $sampleID
echo "rm $CDT3_sense" >> $sampleID
echo "rm $CDT3_anti" >> $sampleID
echo "rm $OUT4" >> $sampleID
echo "rm $CDT4" >> $sampleID
echo "rm $CDT4_sense_gz" >> $sampleID
echo "rm $CDT4_anti_gz" >> $sampleID
echo "rm $CDT4_sense" >> $sampleID
echo "rm $CDT4_anti" >> $sampleID
echo "rm $OUT5" >> $sampleID
echo "rm $CDT5" >> $sampleID
echo "rm $CDT5_sense_gz" >> $sampleID
echo "rm $CDT5_anti_gz" >> $sampleID
echo "rm $CDT5_sense" >> $sampleID
echo "rm $CDT5_anti" >> $sampleID
echo "rm $NT_count" >> $sampleID
echo "rm $MASKED_region" >> $sampleID
echo "rm $category1_sense_smoothed_3" >> $sampleID
echo "rm $category1_anti_smoothed_3" >> $sampleID
echo "rm $category2_sense_smoothed_3" >> $sampleID
echo "rm $category2_anti_smoothed_3" >> $sampleID
echo "rm $category3_sense_smoothed_3" >> $sampleID
echo "rm $category3_anti_smoothed_3" >> $sampleID
echo "rm $category4_sense_smoothed_3" >> $sampleID
echo "rm $category4_anti_smoothed_3" >> $sampleID
echo "rm $category1_sense_smoothed_20" >> $sampleID
echo "rm $category1_anti_smoothed_20" >> $sampleID
echo "rm $category2_sense_smoothed_20" >> $sampleID
echo "rm $category2_anti_smoothed_20" >> $sampleID
echo "rm $category3_sense_smoothed_20" >> $sampleID
echo "rm $category3_anti_smoothed_20" >> $sampleID
echo "rm $category4_sense_smoothed_20" >> $sampleID
echo "rm $category4_anti_smoothed_20" >> $sampleID
echo "rm $category1_sense_max" >> $sampleID
echo "rm $category2_sense_max" >> $sampleID
echo "rm $category3_sense_max" >> $sampleID
echo "rm $category4_sense_max" >> $sampleID
echo "rm $all_max_values" >> $sampleID
echo "rm $scale_values" >> $sampleID
echo "rm $translational_category1_sense" >> $sampleID
echo "rm $translational_category2_sense" >> $sampleID
echo "rm $translational_category3_sense" >> $sampleID
echo "rm $translational_category4_sense" >> $sampleID
echo "rm $translational_category1_anti" >> $sampleID
echo "rm $translational_category2_anti" >> $sampleID
echo "rm $translational_category3_anti" >> $sampleID
echo "rm $translational_category4_anti" >> $sampleID
echo "rm $translational_category1" >> $sampleID
echo "rm $translational_category2" >> $sampleID
echo "rm $translational_category3" >> $sampleID
echo "rm $translational_category4" >> $sampleID
echo "rm $translational_values" >> $sampleID
echo "rm $PERIOD" >> $sampleID
echo "rm $rotational_category1_sense" >> $sampleID
echo "rm $rotational_category1_anti" >> $sampleID
echo "rm $rotational_category1_sense2" >> $sampleID
echo "rm $rotational_category1_anti2" >> $sampleID
echo "rm $rotational_category1_magnitude" >> $sampleID
echo "rm $rotational_values" >> $sampleID
echo "rm $TSV_ratio" >> $sampleID
echo "#rm q1_nucleosome_region_sense.tab" >> $sampleID
echo "#rm q2_nucleosome_region_sense.tab" >> $sampleID
echo "#rm q3_nucleosome_region_sense.tab" >> $sampleID
echo "#rm q4_nucleosome_region_sense.tab" >> $sampleID
echo "#rm significant_peaks_sense.tab" >> $sampleID
echo "#rm significant_peaks_sense_mode.tab" >> $sampleID
echo "#rm output_filtered_nucleosome_sense.tab" >> $sampleID
echo "#rm shifted_columns_sense.tab" >> $sampleID
echo "#rm category1_sense_smoothed_3_full.tab" >> $sampleID
echo "#rm category2_sense_smoothed_3_full.tab" >> $sampleID
echo "#rm category3_sense_smoothed_3_full.tab" >> $sampleID
echo "#rm category4_sense_smoothed_3_full.tab" >> $sampleID
echo "#rm q1_nucleosome_region_anti.tab" >> $sampleID
echo "#rm q2_nucleosome_region_anti.tab" >> $sampleID
echo "#rm q3_nucleosome_region_anti.tab" >> $sampleID
echo "#rm q4_nucleosome_region_anti.tab" >> $sampleID
echo "#rm significant_peaks_anti.tab" >> $sampleID
echo "#rm significant_peaks_anti_mode.tab" >> $sampleID
echo "#rm output_filtered_nucleosome_anti.tab" >> $sampleID
echo "#rm shifted_columns_anti.tab" >> $sampleID
echo "#rm category1_anti_smoothed_3_full.tab" >> $sampleID
echo "#rm category2_anti_smoothed_3_full.tab" >> $sampleID
echo "#rm category3_anti_smoothed_3_full.tab" >> $sampleID
echo "#rm category4_anti_smoothed_3_full.tab" >> $sampleID
echo "#rm category1_rotational_magnitude.tab" >> $sampleID
echo "#rm category2_rotational_magnitude.tab" >> $sampleID
echo "#rm category3_rotational_magnitude.tab" >> $sampleID
echo "#rm category4_rotational_magnitude.tab" >> $sampleID
echo "#rm SENSE_count.tab" >> $sampleID
echo "#rm ANTI_count.tab" >> $sampleID
echo "#script DONE" >> $sampleID
