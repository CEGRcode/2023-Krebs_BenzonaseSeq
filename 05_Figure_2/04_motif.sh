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
ENCODE_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/files/ENCFF868QLL.bed.gz
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/final_bedfiles/MA0834_1_final_1000bp.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/231017_Encode_motif/04_DNAshape_output

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
ENCODE_BEDFILE_unzipped=$(echo $ENCODE_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1".bed"}')
ENCODE_BEDFILE_shuffled=$(echo $ENCODE_BEDFILE_unzipped | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
BEDFILE_shuffled=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
ENCODE_BEDFILE_1000bp=$(echo $ENCODE_BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BEDFILE_20bp=$(echo $BEDFILE_shuffled | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_20bp.bed"}')
TARGET_INTERSECT=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_intersected.bed"}')
NUMBER=$(echo "$TARGET_INTERSECT" | awk -F. '{print $1"_rowsNumber.tab"}')
TARGET_lowlyBound=$(echo "$TARGET_INTERSECT" | awk -F. '{print $1"_lowlyBound.bed"}')
TARGET_lowlyBound_1000bp=$(echo $TARGET_lowlyBound | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_1000bp.bed"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads.out"}')
CDT1=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads"}')
CDT1_sense_gz=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads_sense.cdt.gz"}')
CDT1_anti_gz=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads_anti.cdt.gz"}')
CDT1_sense=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads_sense.cdt"}')
CDT1_anti=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads_anti.cdt"}')
SCALE1=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_ForCDT_allReads"}')
SCALE1a=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_ForCDT_allReads_ScalingFactors.out"}')
CDT1_SCALED_sense=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads_sense_scaled.cdt"}')
CDT1_SCALED_anti=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_allReads_anti_scaled.cdt"}')
SCALED_OUT1_sense=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_sense.tab"}')
SCALED_OUT1_anti=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_ForComposite_scaled_allReads_anti.tab"}')
SCALED_OUT1_final=$(echo "$BAM1a""_""$TARGET_lowlyBound_1000bp" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')
DNAshape_lowlyBound=$(echo $TARGET_lowlyBound_1000bp | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')

sampleID=ATF7_DNAshape_lowlyBound_v4_231202.slurm
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
echo "bedtools intersect -wo -a $ENCODE_BEDFILE_1000bp -b $BEDFILE_20bp -bed > $TARGET_INTERSECT" >> $sampleID
echo "#get number of rows from intersected bedfile" >> $sampleID
echo "cat $TARGET_INTERSECT | wc -l | awk '{printf \"%.f\\n\", \$1 * 0.5}' > $NUMBER" >> $sampleID
echo "#sort intersected bedfile by ChIP enrichment value ("signal value") column 7 (from Encode peaks bedfile), then make bedfile from columns 11-16 (from motif bedfile)." >> $sampleID
echo "cat $TARGET_INTERSECT | sed 's/ //g' | sort -k7,7rn | tail -\$(cat $NUMBER ) | awk '{print \$11\"\t\"\$12\"\t\"\$13\"\t\"\$14\"\t\"\$15\"\t\"\$16}' > $TARGET_lowlyBound" >> $sampleID
echo "#expand befiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=1000 $TARGET_lowlyBound -o=$TARGET_lowlyBound_1000bp" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -a -5 -z --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $TARGET_lowlyBound_1000bp $BAM1" >> $sampleID
echo "#scale output files: options - total tag scaling -t" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1 $BAM1" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1_sense_gz > $CDT1_sense" >> $sampleID
echo "gunzip -c $CDT1_anti_gz > $CDT1_anti" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_sense --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_anti --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_anti" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT1_SCALED_sense $SCALED_OUT1_sense" >> $sampleID
echo "perl $JOB $CDT1_SCALED_anti $SCALED_OUT1_anti" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library." >> $sampleID
echo "cat $SCALED_OUT1_sense $SCALED_OUT1_anti | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT1_final" >> $sampleID
echo "#determine DNA shape with scriptmanager" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER sequence-analysis dna-shape-bed --all --avg-composite -o=$DNAshape_lowlyBound $GENOME $TARGET_lowlyBound_1000bp" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $ENCODE_BEDFILE_unzipped" >> $sampleID
echo "rm $ENCODE_BEDFILE_1000bp" >> $sampleID
echo "rm $BEDFILE_20bp" >> $sampleID
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1" >> $sampleID
echo "rm $CDT1_sense_gz" >> $sampleID
echo "rm $CDT1_anti_gz" >> $sampleID
echo "rm $CDT1_sense" >> $sampleID
echo "rm $CDT1_anti" >> $sampleID
echo "rm $SCALE1a" >> $sampleID
echo "rm $CDT1_SCALED_sense" >> $sampleID
echo "rm $CDT1_SCALED_anti" >> $sampleID
echo "rm $SCALED_OUT1_sense" >> $sampleID
echo "rm $SCALED_OUT1_anti" >> $sampleID
echo "#script DONE" >> $sampleID
