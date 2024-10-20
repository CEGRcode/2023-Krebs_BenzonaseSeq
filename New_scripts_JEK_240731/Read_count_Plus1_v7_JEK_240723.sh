# purpose - Count tags in benz-seq OR MNase-seq at +1 positions / divide by avg. read count / nuc (15.103) -> violin plot out BNase vs. MNase counts at CpG vs non-CpG +1 positions; v3 - scales read counts by seq. depth; v4 - use 5' read1 for counts; v5 - use only +1 positions that are from full-length benz peaks; v6 - math changed for updated MNase seq peaks; v7 - log tranform data before putting into Seaborn, now 240723

# usage
# qq
#
# example
#
# 'qq'

#set bedfiles
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/more_subnucleosomes/fig4c_particles/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet_FullLengthNuc_Only.bed
CpG_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/bedfiles/UCSCgb_hg19_CpGislands_230426.bed

#Scaling files
SCALE1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240702_Quantify/FILES/K562_benzonase-seq_master_midpoint_ForCDT_ScalingFactors.out
SCALE2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240702_Quantify/FILES/K562_MNase_midpoint_ForCDT_ScalingFactors.out

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240702_Quantify/OUTPUT_v7

#set bam files
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/K562_benzonase-seq_master.bam
BAM2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/MNase/K562_MNase.bam

#set genome and human.hg19.genome file
GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/referenceDATA_Will/GENOMES/hg19.fa
HG19_GENOME=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_master_bedfile/files/human.hg19.genome

#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl
PLOT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/240702_Quantify/JOB/violin_plots_mod3_v3_240703.py

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

module load anaconda #configures shell to use conda activate
conda activate plot"

#set output file names
BEDFILE_shuffled=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_shuffled.bed"}')
TARGET_INTERSECT=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_intersected.bed"}')
TARGET_noINTERSECT=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_NOintersect.bed"}')
NUMBER_INTERSECT=$(echo "$TARGET_INTERSECT" | awk -F. '{print $1"_rowsNumber.tab"}')
NUMBER_noINTERSECT=$(echo "$TARGET_noINTERSECT" | awk -F. '{print $1"_rowsNumber.tab"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM1b=$(echo $BAM1a | rev | cut -d"_" -f2 | rev | awk -F. '{print $1}')
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BAM2b=$(echo $BAM2a | rev | cut -d"_" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_read1.out"}')
CDT1=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1}')
CDT1b=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT1c=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined.cdt"}')
OUT2=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_read1.out"}')
CDT2=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1}')
CDT2b=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT2c=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined.cdt"}')
OUT3=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_read1.out"}')
CDT3=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1}')
CDT3b=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT3c=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined.cdt"}')
OUT4=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_read1.out"}')
CDT4=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1}')
CDT4b=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT4c=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined.cdt"}')
SCALE1a=$(echo "$BAM1a" | awk -F. '{print $1"_midpoint_ForCDT_ScalingFactors.out"}')
SCALE2a=$(echo "$BAM2a" | awk -F. '{print $1"_midpoint_ForCDT_ScalingFactors.out"}')
CDT1_SCALED=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_scaled.cdt"}')
CDT2_SCALED=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_scaled.cdt"}')
CDT3_SCALED=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_scaled.cdt"}')
CDT4_SCALED=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_scaled.cdt"}')
CDT1_combined_sum=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined_sum.tsv"}')
CDT2_combined_sum=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined_sum.tsv"}')
CDT3_combined_sum=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined_sum.tsv"}')
CDT4_combined_sum=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined_sum.tsv"}')
CDT1_combined_sum_noHeader=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined_sum_noHeader.tsv"}')
CDT2_combined_sum_noHeader=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_combined_sum_noHeader.tsv"}')
CDT3_combined_sum_noHeader=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined_sum_noHeader.tsv"}')
CDT4_combined_sum_noHeader=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_combined_sum_noHeader.tsv"}')
CDT1_read_avg=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_read_avg.tsv"}')
CDT2_read_avg=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_read_avg.tsv"}')
CDT3_read_avg=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_read_avg.tsv"}')
CDT4_read_avg=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_read_avg.tsv"}')
CDT1_read_avg_ID=$(echo "$BAM1a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_read_avg_ID.tsv"}')
CDT2_read_avg_ID=$(echo "$BAM2a""_""$NUMBER_INTERSECT" | awk -F. '{print $1"_read_avg_ID.tsv"}')
CDT3_read_avg_ID=$(echo "$BAM1a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_read_avg_ID.tsv"}')
CDT4_read_avg_ID=$(echo "$BAM2a""_""$NUMBER_noINTERSECT" | awk -F. '{print $1"_read_avg_ID.tsv"}')
final_readCounts_input=$(echo "readCounts" | awk -F. '{print $1"_file.tsv"}')
SVG=$(echo "readCounts" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1".svg"}')

sampleID=Read_counts_plus1_ROAR_FLNonly_v7_240723.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#shuffle bedfiles" >> $sampleID
echo "shuf $BEDFILE > $BEDFILE_shuffled" >> $sampleID
echo "#do NOT expand ANY befiles" >> $sampleID
echo "#intersect; wo -> so file can be sorted by be intersected if necessary" >> $sampleID
echo "bedtools intersect -u -a $BEDFILE_shuffled -b $CpG_BEDFILE -bed > $TARGET_INTERSECT" >> $sampleID
echo "bedtools intersect -v -a $BEDFILE_shuffled -b $CpG_BEDFILE -bed > $TARGET_noINTERSECT" >> $sampleID
echo "#get number of rows from intersected bedfile" >> $sampleID
echo "cat $TARGET_INTERSECT | wc -l | awk '{printf \"%.f\\n\", \$1}' > $NUMBER_INTERSECT" >> $sampleID
echo "cat $TARGET_noINTERSECT | wc -l | awk '{printf \"%.f\\n\", \$1}' > $NUMBER_noINTERSECT" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --combined --output-matrix=$CDT1 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $TARGET_INTERSECT $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --combined --output-matrix=$CDT2 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $TARGET_INTERSECT $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --combined --output-matrix=$CDT3 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $TARGET_noINTERSECT $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --combined --output-matrix=$CDT4 -N --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 $TARGET_noINTERSECT $BAM2" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1b > $CDT1c" >> $sampleID
echo "gunzip -c $CDT2b > $CDT2c" >> $sampleID
echo "gunzip -c $CDT3b > $CDT3c" >> $sampleID
echo "gunzip -c $CDT4b > $CDT4c" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT1c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT2c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED --scaling-factor=\$(cat $SCALE1 | cut -f2 | tail -1 | awk '{print \$1}') $CDT3c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED --scaling-factor=\$(cat $SCALE2 | cut -f2 | tail -1 | awk '{print \$1}') $CDT4c" >> $sampleID
echo "#sum the number of tags by each row" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT1_combined_sum -r=1 $CDT1_SCALED" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT2_combined_sum -r=1 $CDT2_SCALED" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT3_combined_sum -r=1 $CDT3_SCALED" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis aggregate-data --sum -m -l=3 -o=$CDT4_combined_sum -r=1 $CDT4_SCALED" >> $sampleID
echo "#remove header line from above file" >> $sampleID
echo "paste $CDT1_combined_sum | awk 'NR>1' > $CDT1_combined_sum_noHeader" >> $sampleID
echo "paste $CDT2_combined_sum | awk 'NR>1' > $CDT2_combined_sum_noHeader" >> $sampleID
echo "paste $CDT3_combined_sum | awk 'NR>1' > $CDT3_combined_sum_noHeader" >> $sampleID
echo "paste $CDT4_combined_sum | awk 'NR>1' > $CDT4_combined_sum_noHeader" >> $sampleID
echo "#calculate ratio - updated and checked code: awk 'BEGIN {OFS=\"\t\"}{x = \$2; z = log((x+1)/(15.104+1))/log(10); printf \"%.3f\n\", z}'; Benzonase number is for all nucleosomes" >> $sampleID
echo "cat $CDT1_combined_sum_noHeader | awk 'BEGIN {OFS=\"\t\"}{x = \$2; z = log((x+1)/(15.104+1))/log(10); printf \"%.3f\n\", z}' > $CDT1_read_avg" >> $sampleID
echo "cat $CDT2_combined_sum_noHeader | awk 'BEGIN {OFS=\"\t\"}{x = \$2; z = log((x+1)/(89.660+1))/log(10); printf \"%.3f\n\", z}' > $CDT2_read_avg" >> $sampleID
echo "cat $CDT3_combined_sum_noHeader | awk 'BEGIN {OFS=\"\t\"}{x = \$2; z = log((x+1)/(15.104+1))/log(10); printf \"%.3f\n\", z}' > $CDT3_read_avg" >> $sampleID
echo "cat $CDT4_combined_sum_noHeader | awk 'BEGIN {OFS=\"\t\"}{x = \$2; z = log((x+1)/(89.660+1))/log(10); printf \"%.3f\n\", z}' > $CDT4_read_avg" >> $sampleID
echo "#add column with ID of file" >> $sampleID
echo "cat $CDT1_read_avg | awk '{print \$1\"\t\"\"$BAM1b\"\"_CpG\"}' > $CDT1_read_avg_ID" >> $sampleID
echo "cat $CDT2_read_avg | awk '{print \$1\"\t\"\"$BAM2b\"\"_CpG\"}' > $CDT2_read_avg_ID" >> $sampleID
echo "cat $CDT3_read_avg | awk '{print \$1\"\t\"\"$BAM1b\"\"_non-CpG\"}' > $CDT3_read_avg_ID" >> $sampleID
echo "cat $CDT4_read_avg | awk '{print \$1\"\t\"\"$BAM2b\"\"_non-CpG\"}' > $CDT4_read_avg_ID" >> $sampleID
echo "#make file of all data" >> $sampleID
echo "cat $CDT1_read_avg_ID $CDT2_read_avg_ID $CDT3_read_avg_ID $CDT4_read_avg_ID > $final_readCounts_input" >> $sampleID
echo "#make violin plot wity modified violin_plots_mod.py" >> $sampleID
echo "python $PLOT -i $final_readCounts_input -o $SVG" >> $sampleID
echo "# finish script" >> $sampleID
