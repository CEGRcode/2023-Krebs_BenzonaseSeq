# purpose - intersect bedfile of all sites for a motif with encode-called peaks (Bed narrowPeak bedfile). THEN designate sites that do and do NOT overlap with our 'true NFR' bedfile. **Update 231128 with updated 'TRUE NFR' bedfile

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
ENCODE_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/files/ENCFF868QLL.bed.gz
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/final_bedfiles/MA0834_1_final_1000bp.bed
NFR_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/02_NFR_output_231128/K562_trueNFR.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231128_trueNFR_motif/01_trueNFR_output_231128

#set bam library file to BI_rep1 and set ATACseq library
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/K562_benzonase-seq_master.bam

#set blacklist and .genome file
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed
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
#SBATCH --time=4:00:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

#set output file names
ENCODE_BEDFILE_unzipped=$(echo $ENCODE_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1".bed"}')
ENCODE_BEDFILE_shuffled=$(echo $ENCODE_BEDFILE_unzipped | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
BEDFILE_shuffled=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
ENCODE_BEDFILE_1000bp=$(echo $ENCODE_BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BEDFILE_20bp=$(echo $BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_20bp.bed"}')
TARGET_BOUND=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_bound.bed"}')
NUMBER_BOUND=$(echo "$TARGET_BOUND" | awk -F. '{print $1"_rowsNumber.tab"}')
TARGET_NFR=$(echo $TARGET_BOUND | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NFR.bed"}')
TARGET_nonNFR=$(echo $TARGET_BOUND | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_nonNFR.bed"}')
NUMBER_BOUND_NFR=$(echo "$TARGET_NFR" | awk -F. '{print $1"_rowsNumber.tab"}')
NUMBER_BOUND_nonNFR=$(echo "$TARGET_nonNFR" | awk -F. '{print $1"_rowsNumber.tab"}')
NUMBER_OUTPUT=$(echo "$BEDFILE" | awk -F. '{print $1"_output.tab"}')
TARGET_BOUND_NFR_1000bp=$(echo $TARGET_NFR | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
TARGET_BOUND_nonNFR_1000bp=$(echo $TARGET_nonNFR | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT1=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads"}')
CDT1_sense_gz=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT1_anti_gz=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT1_sense=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT1_anti=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
SCALE1=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_ForCDT_allReads"}')
SCALE1a=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_ForCDT_allReads_ScalingFactors.out"}')
CDT1_SCALED_sense=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads_sense_scaled.cdt"}')
CDT1_SCALED_anti=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_allReads_anti_scaled.cdt"}')
SCALED_OUT1_sense=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_sense.tab"}')
SCALED_OUT1_anti=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_anti.tab"}')
SCALED_OUT1_final=$(echo "$BAM1a""_""$TARGET_BOUND_NFR_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
OUT2=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT2=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads"}')
CDT2_sense_gz=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT2_anti_gz=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT2_sense=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT2_anti=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
SCALE2=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_ForCDT_allReads"}')
SCALE2a=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_ForCDT_allReads_ScalingFactors.out"}')
CDT2_SCALED_sense=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads_sense_scaled.cdt"}')
CDT2_SCALED_anti=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_allReads_anti_scaled.cdt"}')
SCALED_OUT2_sense=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_sense.tab"}')
SCALED_OUT2_anti=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_anti.tab"}')
SCALED_OUT2_final=$(echo "$BAM1a""_""$TARGET_BOUND_nonNFR_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')


sampleID=ATF7_01_trueNFRs_PileUp_v2_231128.slurm
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
echo "#intersect" >> $sampleID
echo "bedtools intersect -wb -a $ENCODE_BEDFILE_1000bp -b $BEDFILE_20bp -bed > $TARGET_BOUND" >> $sampleID
echo "#get number of rows from intersected bedfile" >> $sampleID
echo "cat $TARGET_BOUND | wc -l | awk '{printf \"%.f\\n\", \$1}' > $NUMBER_BOUND" >> $sampleID
echo "#intersect bedfile of bound target with ATACseq reads. These are bound sites in NFRs.*these currently have a 20 bp width." >> $sampleID
echo "bedtools intersect -u -a $TARGET_BOUND -b $NFR_BEDFILE -bed > $TARGET_NFR" >> $sampleID
echo "#determine what bound sites do not overlap with ATACseq reads. These are bound and NOT in NFRs." >> $sampleID
echo "bedtools intersect -v -a $TARGET_BOUND -b $NFR_BEDFILE -bed > $TARGET_nonNFR" >> $sampleID
echo "#get number of rows from NFR and nonNFR bedfile" >> $sampleID
echo "cat $TARGET_NFR | wc -l | awk '{printf \"%.f\\n\", \$1}' > $NUMBER_BOUND_NFR" >> $sampleID
echo "cat $TARGET_nonNFR | wc -l | awk '{printf \"%.f\\n\", \$1}' > $NUMBER_BOUND_nonNFR" >> $sampleID
echo "#expand befiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $TARGET_NFR -o=$TARGET_BOUND_NFR_1000bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $TARGET_nonNFR -o=$TARGET_BOUND_nonNFR_1000bp" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $TARGET_BOUND_NFR_1000bp $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $TARGET_BOUND_nonNFR_1000bp $BAM1" >> $sampleID
echo "#scale output files: options - total tag scaling -t" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2 $BAM1" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1_sense_gz > $CDT1_sense" >> $sampleID
echo "gunzip -c $CDT1_anti_gz > $CDT1_anti" >> $sampleID
echo "gunzip -c $CDT2_sense_gz > $CDT2_sense" >> $sampleID
echo "gunzip -c $CDT2_anti_gz > $CDT2_anti" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_sense --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_anti --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_sense --scaling-factor=\$(cat $SCALE2a | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_anti --scaling-factor=\$(cat $SCALE2a | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_anti" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT1_SCALED_sense $SCALED_OUT1_sense" >> $sampleID
echo "perl $JOB $CDT1_SCALED_anti $SCALED_OUT1_anti" >> $sampleID
echo "perl $JOB $CDT2_SCALED_sense $SCALED_OUT2_sense" >> $sampleID
echo "perl $JOB $CDT2_SCALED_anti $SCALED_OUT2_anti" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library." >> $sampleID
echo "cat $SCALED_OUT1_sense $SCALED_OUT1_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT1_final" >> $sampleID
echo "cat $SCALED_OUT2_sense $SCALED_OUT2_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT2_final" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $ENCODE_BEDFILE_unzipped" >> $sampleID
echo "rm $ENCODE_BEDFILE_shuffled" >> $sampleID
echo "rm $BEDFILE_shuffled" >> $sampleID
echo "rm $ENCODE_BEDFILE_1000bp" >> $sampleID
echo "rm $BEDFILE_20bp" >> $sampleID
echo "rm $NUMBER_OUTPUT" >> $sampleID
echo "rm $TARGET_BOUND_NFR_1000bp" >> $sampleID
echo "rm $TARGET_BOUND_nonNFR_1000bp" >> $sampleID
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1" >> $sampleID
echo "rm $CDT1_sense_gz" >> $sampleID
echo "rm $CDT1_anti_gz" >> $sampleID
echo "rm $CDT1_sense" >> $sampleID
echo "rm $CDT1_anti" >> $sampleID
echo "rm $SCALE1" >> $sampleID
echo "rm $SCALE1a" >> $sampleID
echo "rm $CDT1_SCALED_sense" >> $sampleID
echo "rm $CDT1_SCALED_anti" >> $sampleID
echo "rm $SCALED_OUT1_sense" >> $sampleID
echo "rm $SCALED_OUT1_anti" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2" >> $sampleID
echo "rm $CDT2_sense_gz" >> $sampleID
echo "rm $CDT2_anti_gz" >> $sampleID
echo "rm $CDT2_sense" >> $sampleID
echo "rm $CDT2_anti" >> $sampleID
echo "rm $SCALE2" >> $sampleID
echo "rm $SCALE2a" >> $sampleID
echo "rm $CDT2_SCALED_sense" >> $sampleID
echo "rm $CDT2_SCALED_anti" >> $sampleID
echo "rm $SCALED_OUT2_sense" >> $sampleID
echo "rm $SCALED_OUT2_anti" >> $sampleID
echo "#script DONE" >> $sampleID
