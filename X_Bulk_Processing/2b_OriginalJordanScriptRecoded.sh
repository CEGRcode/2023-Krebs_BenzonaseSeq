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
MEMEFILE=../data/JASPAR/${TARGET}_${JASPAR}.meme
BAMFILE=../data/BAM/BNase-seq_50U-10min_merge_hg38.bam		#BAM1
BLACKLIST=../data/hg38_files/ENCFF356LFX_hg38_exclude.bed.gz

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240909_TFBS/01_final_TF_pipeline_jobs_240913/CTCF_NucOccupancy_settings_pipeline_MA1929_1_v13_241011

# Script shortcuts
SCRIPTMANAGER=../bin/ScriptManager-v0.15.jar
DEDUP=../bin/dedup_coord_by_ID.py
EXTRACT=../bin/extract_row_number_240817.py
MASKED=../bin/masked_region_240817.py
SMOOTH=../bin/smoothing_parameterize.py
MAX=../bin/max_position_v3_240818.py
SCALE=../bin/scaling_240814.py
TRANSLATIONAL=../bin/translational_range_parameterize.py
TRANSLATIONAL_average=../bin/translational_range_average_240820.py
AUTO=../bin/autocorrelation_of_CDT_v2_240818.py
PERIODICITY=../bin/periodicity_240818.py
ROTATIONAL=../bin/rotational_ratio_parameterize.py
MODE_sense=../bin/rotational_sense_mode_v2_240826.py
MODE_sense_substitute=../bin/MODE_sense_substitute_241011.py
PEAKS_shift=../bin/rotational_peaks_shift_v2_241011.py
PEAKS_fill=../bin/rotational_peaks_shifted_columns_v3_240825.py
MODE_anti=../bin/rotational_anti_mode_v2_240826.py
FILTER=../bin/rotational_peaks_filter_240826.py
ROTATIONAL_magnitude=../bin/rotational_magnitude_v2_240826.py
SENSE_count=../bin/rotational_sense_count_240826.py
ANTI_count=../bin/rotational_anti_count_240826.py
CONCAT=../bin/concatenate_v3_241011.py

# ===============================================================================================================================

RUNID=${TARGET}_${JASPAR}
MOTIF=../data/RefPT-JASPAR

#set output file names
BEDFILE_category1=${RUNID}_SORT-TFnucRatio_GROUP-Quartile1
BEDFILE_category2=${RUNID}_SORT-TFnucRatio_GROUP-Quartile2
BEDFILE_category3=${RUNID}_SORT-TFnucRatio_GROUP-Quartile3
BEDFILE_category4=${RUNID}_SORT-TFnucRatio_GROUP-Quartile4
BEDFILE_category1_1000bp=$MOTIF/1000bp/${RUNID}_SORT-TFnucRatio_GROUP-Quartile1_1000bp.bed
BEDFILE_category2_1000bp=$MOTIF/1000bp/${RUNID}_SORT-TFnucRatio_GROUP-Quartile2_1000bp.bed
BEDFILE_category3_1000bp=$MOTIF/1000bp/${RUNID}_SORT-TFnucRatio_GROUP-Quartile3_1000bp.bed
BEDFILE_category4_1000bp=$MOTIF/1000bp/${RUNID}_SORT-TFnucRatio_GROUP-Quartile4_1000bp.bed

NT_count=${RUNID}_NT_count.tab
MASKED_region=${RUNID}_masked.tab

category1_sense_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_ForComposite_allReads_sense_smooth3.tab
category2_sense_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_ForComposite_allReads_sense_smooth3.tab
category3_sense_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_ForComposite_allReads_sense_smooth3.tab
category4_sense_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_ForComposite_allReads_sense_smooth3.tab
category1_anti_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_ForComposite_allReads_anti_smooth3.tab
category2_anti_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_ForComposite_allReads_anti_smooth3.tab
category3_anti_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_ForComposite_allReads_anti_smooth3.tab
category4_anti_smoothed_3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_ForComposite_allReads_anti_smooth3.tab

scale_values=${RUNID}_scale_values.tab

translational_category1=${RUNID}_GROUP-Quartile1_translational_setting.tab
translational_category2=${RUNID}_GROUP-Quartile2_translational_setting.tab
translational_category3=${RUNID}_GROUP-Quartile3_translational_setting.tab
translational_category4=${RUNID}_GROUP-Quartile4_translational_setting.tab

translational_values=${RUNID}_translational_values.tab
PERIOD=${RUNID}_periodicity.tab
correlation_results=${RUNID}_correlation_results


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

# See 03_Call_Motifs for generating initial motifs split into quartiles

#extract number of NTs from MEME file
python $EXTRACT $MEMEFILE $NT_count
#determine the 5' and 3' boundaries of the motif masked region relative to the center column of tab files at column 501
python $MASKED $NT_count $MASKED_region

for QUARTILE in {1..4};
do
	CATEGORY=${RUNID}_SORT-TFnucRatio_GROUP-Quartile${QUARTILE}
	BEDFILE=$MOTIF/1000bp/${RUNID}_SORT-TFnucRatio_GROUP-Quartile${QUARTILE}_1000bp.bed
	OUT_COMPOSITE=0${QUARTILE}_BNase-seq_50U-10min_merge_hg38_${CATEGORY}_ForComposite_final.tab

	BASE=BNase-seq_50U-10min_merge_hg38_${CATEGORY}_ForComposite_allReads

	# Tag pileup Benzonase cut-sites. Settings: 5 prime end (-5) with read 1 (-1), No smoothing (N), required proper PEs (-p)
	java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -N -p --cpu 4 -o $OUT_COMPOSITE $BEDFILE $BAMFILE

	# Slice sense and anti strand
	awk 'NR==1;NR==2' $OUT_COMPOSITE > ${BASE}_sense.tab
	awk 'NR==1;NR==3' $OUT_COMPOSITE > ${BASE}_anti.tab

	# Apply smoothing (3) to sense and anti
	python $SMOOTH 3 ${BASE}_sense.tab ${BASE}_sense_smooth3.tab	#category1_sense_smoothed_3
	python $SMOOTH 3 ${BASE}_anti.tab ${BASE}_anti_smooth3.tab		#category1_anti_smoothed_3

	# Apply smoothing (20) to sense and anti
	python $SMOOTH 20 ${BASE}_sense.tab ${BASE}_sense_smooth20.tab
	python $SMOOTH 20 ${BASE}_anti.tab ${BASE}_anti_smooth20.tab

	# Get max positions (for later scaling) of sense strand from column 276 (bp -225) - 326 (bp-175) AND determine the bp of the max position. OUTPUT file is name, max value, position of max value
	python $MAX ${BASE}_sense_smooth20.tab ${BASE}_sense_smooth20_max.tab

	# === Check translational setting ===

	RBASE=${RUNID}_GROUP-Quartile${QUARTILE}

	#get max range (max-min) from -350 to -150 bp for motif strand
	python $TRANSLATIONAL sense ${BASE}_sense_smooth20.tab ${RBASE}_sense_translational_setting.tab

	#get max range (max-min) from +150 to +350 bp for opposite strand
	python $TRANSLATIONAL anti ${BASE}_sense_smooth20.tab ${RBASE}_anti_translational_setting.tab

	#get average of range of translational magnitude from both strands
	python $TRANSLATIONAL_average ${RBASE}_sense_translational_setting.tab ${RBASE}_anti_translational_setting.tab ${RBASE}_translational_setting.tab #translational_category1

	# === Check rotation ===

	#get significant peaks from category1 sense strand and use those respective bins to call peaks from sense strands of categories 2, 3, and 4
	python $ROTATIONAL sense ${BASE}_sense_smooth3.tab ${RUNID}_q${QUARTILE}_nucleosome_region_sense.tab

	#get significant peaks from category1 anti strand and use those respective bins to call peaks from sense strands of categories 2, 3, and 4
	python $ROTATIONAL anti ${BASE}_sense_smooth3.tab ${RUNID}_q${QUARTILE}_nucleosome_region_anti.tab

done


Q1=BNase-seq_50U-10min_merge_hg38_${RUNID}_SORT-TFnucRatio_GROUP-Quartile1_ForComposite_allReads
Q2=BNase-seq_50U-10min_merge_hg38_${RUNID}_SORT-TFnucRatio_GROUP-Quartile2_ForComposite_allReads
Q3=BNase-seq_50U-10min_merge_hg38_${RUNID}_SORT-TFnucRatio_GROUP-Quartile3_ForComposite_allReads
Q4=BNase-seq_50U-10min_merge_hg38_${RUNID}_SORT-TFnucRatio_GROUP-Quartile4_ForComposite_allReads

#combine all above tab files (and remove headers of last 3)
cat ${Q1}_sense_smooth20_max.tab \
	${Q2}_sense_smooth20_max.tab \
	${Q3}_sense_smooth20_max.tab \
	${Q4}_sense_smooth20_max.tab \
	| awk 'NR==1;NR==2;NR==4;NR==6;NR==8' \
	> ${RUNID}_all_max_values.tab

#get scaling value for all categories
python $SCALE ${RUNID}_all_max_values.tab $scale_values

#combine all above tab files (and add first column of quartile info and header). Output is average of peaks from either strand
cat ${RUNID}_GROUP-Quartile1_translational_setting.tab \
	${RUNID}_GROUP-Quartile2_translational_setting.tab \
	${RUNID}_GROUP-Quartile3_translational_setting.tab \
	${RUNID}_GROUP-Quartile4_translational_setting.tab \
	| awk 'BEGIN{print "Average_Translational_Magnitude"}1' \
	| awk 'BEGIN{quartile[1]="Quartile"; for(i=2;i<=5;i++) quartile[i]=i-1} {print quartile[NR]"\t"$0} NR>5' \
	> $translational_values

#perform autocorrelation to determine most likely periodicity
python $AUTO -i BNase-seq_50U-10min_merge_hg38_${RUNID}_SORT-TFnucRatio_GROUP-Quartile1_ForComposite_allReads_sense_smooth3.tab -o $correlation_results
python $PERIODICITY $correlation_results.tsv $PERIOD


#get concatenated list of all unique, significant peaks, THEN get their the mode of their max position (bp)
cat ${RUNID}_q1_nucleosome_region_sense.tab \
	${RUNID}_q2_nucleosome_region_sense.tab \
	${RUNID}_q3_nucleosome_region_sense.tab \
	${RUNID}_q4_nucleosome_region_sense.tab \
	> ${RUNID}_significant_peaks_sense.tab
python $MODE_sense ${RUNID}_significant_peaks_sense.tab $MASKED_region ${RUNID}_significant_peaks_sense_mode.tab

#determine unique set of 'significant' peaks from motif strand and do unique sort
cut -f1,2  ${RUNID}_significant_peaks_sense.tab \
	| sort -k1,1 | uniq \
	> ${RUNID}_output_filtered_nucleosome_sense.tab

#Sense mode needs substituted to match 'opposite strand phase (0-9) then shift by doing by 5' or 3' by 'mode-5' with +=shift 5', -=shift3'
python $MODE_sense_substitute ${RUNID}_significant_peaks_sense_mode.tab ${RUNID}_significant_peaks_sense_mode_substituted.tab

#sort unique, significant peaks by the above substituted mode'
python $PEAKS_shift ${RUNID}_output_filtered_nucleosome_sense.tab ${RUNID}_significant_peaks_sense_mode_substituted.tab ${RUNID}_shifted_columns_sense.tab


##repeat above for opposite strand

#get concatenated list of all unique, significant peaks, THEN get their the mode of their max position (bp)
cat ${RUNID}_q1_nucleosome_region_anti.tab \
	${RUNID}_q2_nucleosome_region_anti.tab \
	${RUNID}_q3_nucleosome_region_anti.tab \
	${RUNID}_q4_nucleosome_region_anti.tab \
	> ${RUNID}_significant_peaks_anti.tab
python $MODE_anti ${RUNID}_significant_peaks_anti.tab $MASKED_region ${RUNID}_significant_peaks_anti_mode.tab

#determine unique set of 'significant' peaks from opposite strand and do unique sort
cut -f1,2 ${RUNID}_significant_peaks_anti.tab \
	| sort -k1,1 | uniq \
	> ${RUNID}_output_filtered_nucleosome_anti.tab

#sort unique, significant peaks by the above mode; shift is by adding mode/2 to shift 3' to the motif
python $PEAKS_shift ${RUNID}_output_filtered_nucleosome_anti.tab ${RUNID}_significant_peaks_anti_mode.tab ${RUNID}_shifted_columns_anti.tab



for QUARTILE in {1..4};
do
	BEDFILE=$MOTIF/1000bp/${RUNID}_SORT-TFnucRatio_GROUP-Quartile${QUARTILE}_1000bp.bed
	OUT_COMPOSITE=0${QUARTILE}_BNase-seq_50U-10min_merge_hg38_${CATEGORY}_ForComposite_final.tab

	CATEGORY=${RUNID}_SORT-TFnucRatio_GROUP-Quartile${QUARTILE}
	BASE=BNase-seq_50U-10min_merge_hg38_${CATEGORY}_ForComposite_allReads

	smoothed_base=BNase-seq_50U-10min_merge_hg38_${CATEGORY}_ForComposite_allReads

	#take shifted, significant bins and fill out range for each bin
	python $PEAKS_fill ${smoothed_base}_sense_smooth3.tab ${RUNID}_shifted_columns_sense.tab category${QUARTILE}_sense_smoothed_3_full.tab
	python $PEAKS_fill ${smoothed_base}_anti_smooth3.tab ${RUNID}_shifted_columns_anti.tab category${QUARTILE}_anti_smoothed_3_full.tab

done

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
python $CONCAT $scale_values $PERIOD $translational_values ${RUNID}_significant_peaks_sense_mode_substituted.tab ${RUNID}_significant_peaks_anti_mode.tab $rotational_values SENSE_count.tab ANTI_count.tab $FINAL
