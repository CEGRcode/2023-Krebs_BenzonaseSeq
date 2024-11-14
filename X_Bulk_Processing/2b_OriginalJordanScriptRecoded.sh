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
BAM1a=BNase-seq_50U-10min_merge_hg38

OUT2=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_1000bp_allReads.out
CDT2=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_1000bp_allReads
CDT2_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_allReads_sense.cdt
CDT2_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_allReads_anti.cdt
OUT2_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_ForComposite_allReads_sense.tab
OUT2_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_ForComposite_allReads_anti.tab
OUT2_final=01_BNase-seq_50U-10min_merge_hg38_${BEDFILE_category1}_ForComposite_final.tab

OUT3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_allReads.out
CDT3=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_allReads
CDT3_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_allReads_sense.cdt
CDT3_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_allReads_anti.cdt
OUT3_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_ForComposite_allReads_sense.tab
OUT3_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_ForComposite_allReads_anti.tab
OUT3_final=02_BNase-seq_50U-10min_merge_hg38_${BEDFILE_category2}_ForComposite_final.tab

OUT4=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_allReads.out
CDT4=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_allReads
CDT4_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_allReads_sense.cdt
CDT4_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_allReads_anti.cdt
OUT4_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_ForComposite_allReads_sense.tab
OUT4_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_ForComposite_allReads_anti.tab
OUT4_final=03_BNase-seq_50U-10min_merge_hg38_${BEDFILE_category3}_ForComposite_final.tab

OUT5=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_allReads.out
CDT5=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_allReads
CDT5_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_allReads_sense.cdt
CDT5_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_allReads_anti.cdt
OUT5_sense=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_ForComposite_allReads_sense.tab
OUT5_anti=BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_ForComposite_allReads_anti.tab
OUT5_final=04_BNase-seq_50U-10min_merge_hg38_${BEDFILE_category4}_ForComposite_final.tab

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

# See 03_Call_Motifs for generating initial motifs split into quartiles

#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 --output-matrix=$CDT2 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $BEDFILE_category1_1000bp $BAMFILE
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 --output-matrix=$CDT3 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $BEDFILE_category2_1000bp $BAMFILE
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 --output-matrix=$CDT4 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 $BEDFILE_category3_1000bp $BAMFILE
java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 --output-matrix=$CDT5 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 $BEDFILE_category4_1000bp $BAMFILE
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
