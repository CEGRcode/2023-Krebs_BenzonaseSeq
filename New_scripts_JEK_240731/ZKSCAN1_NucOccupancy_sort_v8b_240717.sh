# purpose - intersect bedfile of all sites for a motif with encode-called peaks (Bed narrowPeak bedfile). v8 - intersection (and not intersection) based on Olivia's code. Rest of code from v6

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
ENCODE_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/files/ENCFF163VUK.bed.gz
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/final_bedfiles/MA1585_1_final_1000bp.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_lowlyBound/ZKSCAN1_NucOccupancy_v8_sort_240717

#set bam library file to DNase-seq SE file; add ENCODE TF BAM file
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/K562_benzonase-seq_master.bam
BAM2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240709_DNase/final_files/SRR14304993_master.bam
BAM3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/MNase/K562_MNase.bam

#Scaling factor files
SCALE1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_lowlyBound/FILES/K562_benzonase-seq_master_ScalingFactors.out
SCALE2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_lowlyBound/FILES/SRR14304993_master_ScalingFactors.out
SCALE3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_lowlyBound/FILES/K562_MNase_ScalingFactors.out

#set blacklist and .genome file
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl
JOB_ROW=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig6_Subnucleosomes/job/sum_Row_CDT.pl
DEDUP=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240703_lowlyBound/jobs/dedup_coord_by_ID.py

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

module load anaconda #configures shell to use conda activate
conda activate bioinfo"

#set output file names
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
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1}')
CDT1b=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT1c=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined.cdt"}')
CDT1_sum=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined_sum.tsv"}')
CDT1_sum_noHeader=$(echo "$BAM1a""_""$TARGET_Bound_164bp" | awk -F. '{print $1"_combined_sum_noHeader.tsv"}')
TSV_ratio=$(echo $TARGET_Bound | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_ratio.TSV"}')
NUMBER=$(echo "$TSV_ratio" | awk -F. '{print $1"_rowsNumber.tab"}')
BEDFILE_category1=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category1.bed"}')
BEDFILE_category2=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category2.bed"}')
BEDFILE_category3=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category3.bed"}')
BEDFILE_category4=$(echo $TARGET_Bound_164bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_category4.bed"}')
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
CDT2_SCALED_sense=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads_sense_scaled.cdt"}')
CDT2_SCALED_anti=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_allReads_anti_scaled.cdt"}')
SCALED_OUT2_sense=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_sense.tab"}')
SCALED_OUT2_anti=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_anti.tab"}')
SCALED_OUT2_final=$(echo "$BAM1a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT3=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT3=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads"}')
CDT3_sense_gz=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT3_anti_gz=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT3_sense=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT3_anti=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
CDT3_SCALED_sense=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_sense_scaled.cdt"}')
CDT3_SCALED_anti=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_allReads_anti_scaled.cdt"}')
SCALED_OUT3_sense=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_sense.tab"}')
SCALED_OUT3_anti=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_anti.tab"}')
SCALED_OUT3_final=$(echo "$BAM1a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT4=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT4=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads"}')
CDT4_sense_gz=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT4_anti_gz=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT4_sense=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT4_anti=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
CDT4_SCALED_sense=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_sense_scaled.cdt"}')
CDT4_SCALED_anti=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_allReads_anti_scaled.cdt"}')
SCALED_OUT4_sense=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_sense.tab"}')
SCALED_OUT4_anti=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_anti.tab"}')
SCALED_OUT4_final=$(echo "$BAM1a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT5=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT5=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads"}')
CDT5_sense_gz=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT5_anti_gz=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT5_sense=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT5_anti=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
CDT5_SCALED_sense=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_sense_scaled.cdt"}')
CDT5_SCALED_anti=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_allReads_anti_scaled.cdt"}')
SCALED_OUT5_sense=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_sense.tab"}')
SCALED_OUT5_anti=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_anti.tab"}')
SCALED_OUT5_final=$(echo "$BAM1a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT6=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1.out"}')
CDT6=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1"}')
CDT6_sense_gz=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1_sense.cdt.gz"}')
CDT6_anti_gz=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1_anti.cdt.gz"}')
CDT6_sense=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1_sense.cdt"}')
CDT6_anti=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1_anti.cdt"}')
CDT6_SCALED_sense=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1_sense_scaled.cdt"}')
CDT6_SCALED_anti=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_Read1_anti_scaled.cdt"}')
SCALED_OUT6_sense=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_sense.tab"}')
SCALED_OUT6_anti=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_anti.tab"}')
SCALED_OUT6_final=$(echo "$BAM2a""_""$BEDFILE_category1_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT7=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1.out"}')
CDT7=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1"}')
CDT7_sense_gz=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1_sense.cdt.gz"}')
CDT7_anti_gz=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1_anti.cdt.gz"}')
CDT7_sense=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1_sense.cdt"}')
CDT7_anti=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1_anti.cdt"}')
CDT7_SCALED_sense=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1_sense_scaled.cdt"}')
CDT7_SCALED_anti=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_Read1_anti_scaled.cdt"}')
SCALED_OUT7_sense=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_sense.tab"}')
SCALED_OUT7_anti=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_anti.tab"}')
SCALED_OUT7_final=$(echo "$BAM2a""_""$BEDFILE_category2_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT8=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1.out"}')
CDT8=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1"}')
CDT8_sense_gz=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1_sense.cdt.gz"}')
CDT8_anti_gz=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1_anti.cdt.gz"}')
CDT8_sense=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1_sense.cdt"}')
CDT8_anti=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1_anti.cdt"}')
CDT8_SCALED_sense=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1_sense_scaled.cdt"}')
CDT8_SCALED_anti=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_Read1_anti_scaled.cdt"}')
SCALED_OUT8_sense=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_sense.tab"}')
SCALED_OUT8_anti=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_anti.tab"}')
SCALED_OUT8_final=$(echo "$BAM2a""_""$BEDFILE_category3_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT9=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1.out"}')
CDT9=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1"}')
CDT9_sense_gz=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1_sense.cdt.gz"}')
CDT9_anti_gz=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1_anti.cdt.gz"}')
CDT9_sense=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1_sense.cdt"}')
CDT9_anti=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1_anti.cdt"}')
CDT9_SCALED_sense=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1_sense_scaled.cdt"}')
CDT9_SCALED_anti=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_Read1_anti_scaled.cdt"}')
SCALED_OUT9_sense=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_sense.tab"}')
SCALED_OUT9_anti=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_scaled_Read1_anti.tab"}')
SCALED_OUT9_final=$(echo "$BAM2a""_""$BEDFILE_category4_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')

sampleID=ZKSCAN1_NucOccupancy_sort_ratio_v8_240717.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
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
echo "cat $TSV_ratio | wc -l | awk '{printf \"%.f\\n\", \$1}' > $NUMBER" >> $sampleID
echo "#make bedfile for sites that have a TF:Nuc ratio of <1, 1-5, 5-10, <10" >> $sampleID
echo "cat $TSV_ratio  | awk '{if (\$7<=1) print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category1" >> $sampleID
echo "cat $TSV_ratio  | awk '{if (\$7>1 && \$7<5) print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category2" >> $sampleID
echo "cat $TSV_ratio  | awk '{if (\$7>5 && \$7<10) print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category3" >> $sampleID
echo "cat $TSV_ratio  | awk '{if (\$7>=10) print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' > $BEDFILE_category4" >> $sampleID
echo "#expand befiles" >> $sampleID
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
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_sense --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_anti --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_sense --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT3_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_anti --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT3_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_sense --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_anti --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_sense --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT5_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_anti --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT5_anti" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT2_SCALED_sense $SCALED_OUT2_sense" >> $sampleID
echo "perl $JOB $CDT2_SCALED_anti $SCALED_OUT2_anti" >> $sampleID
echo "perl $JOB $CDT3_SCALED_sense $SCALED_OUT3_sense" >> $sampleID
echo "perl $JOB $CDT3_SCALED_anti $SCALED_OUT3_anti" >> $sampleID
echo "perl $JOB $CDT4_SCALED_sense $SCALED_OUT4_sense" >> $sampleID
echo "perl $JOB $CDT4_SCALED_anti $SCALED_OUT4_anti" >> $sampleID
echo "perl $JOB $CDT5_SCALED_sense $SCALED_OUT5_sense" >> $sampleID
echo "perl $JOB $CDT5_SCALED_anti $SCALED_OUT5_anti" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library." >> $sampleID
echo "cat $SCALED_OUT2_sense $SCALED_OUT2_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT2_final" >> $sampleID
echo "cat $SCALED_OUT3_sense $SCALED_OUT3_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT3_final" >> $sampleID
echo "cat $SCALED_OUT4_sense $SCALED_OUT4_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT4_final" >> $sampleID
echo "cat $SCALED_OUT5_sense $SCALED_OUT5_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT5_final" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist for DNase-seq**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT6 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6 $BEDFILE_category1_1000bp $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT7 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT7 $BEDFILE_category2_1000bp $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT8 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT8 $BEDFILE_category3_1000bp $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT9 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT9 $BEDFILE_category4_1000bp $BAM2" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT6_sense_gz > $CDT6_sense" >> $sampleID
echo "gunzip -c $CDT6_anti_gz > $CDT6_anti" >> $sampleID
echo "gunzip -c $CDT7_sense_gz > $CDT7_sense" >> $sampleID
echo "gunzip -c $CDT7_anti_gz > $CDT7_anti" >> $sampleID
echo "gunzip -c $CDT8_sense_gz > $CDT8_sense" >> $sampleID
echo "gunzip -c $CDT8_anti_gz > $CDT8_anti" >> $sampleID
echo "gunzip -c $CDT9_sense_gz > $CDT9_sense" >> $sampleID
echo "gunzip -c $CDT9_anti_gz > $CDT9_anti" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_sense --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_anti --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT6_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT7_SCALED_sense --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT7_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT7_SCALED_anti --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT7_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT8_SCALED_sense --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT8_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT8_SCALED_anti --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT8_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT9_SCALED_sense --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT9_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT9_SCALED_anti --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT9_anti" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT6_SCALED_sense $SCALED_OUT6_sense" >> $sampleID
echo "perl $JOB $CDT6_SCALED_anti $SCALED_OUT6_anti" >> $sampleID
echo "perl $JOB $CDT7_SCALED_sense $SCALED_OUT7_sense" >> $sampleID
echo "perl $JOB $CDT7_SCALED_anti $SCALED_OUT7_anti" >> $sampleID
echo "perl $JOB $CDT8_SCALED_sense $SCALED_OUT8_sense" >> $sampleID
echo "perl $JOB $CDT8_SCALED_anti $SCALED_OUT8_anti" >> $sampleID
echo "perl $JOB $CDT9_SCALED_sense $SCALED_OUT9_sense" >> $sampleID
echo "perl $JOB $CDT9_SCALED_anti $SCALED_OUT9_anti" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library." >> $sampleID
echo "cat $SCALED_OUT6_sense $SCALED_OUT6_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT6_final" >> $sampleID
echo "cat $SCALED_OUT7_sense $SCALED_OUT7_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT7_final" >> $sampleID
echo "cat $SCALED_OUT8_sense $SCALED_OUT8_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT8_final" >> $sampleID
echo "cat $SCALED_OUT9_sense $SCALED_OUT9_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT9_final" >> $sampleID
echo "#script DONE" >> $sampleID
