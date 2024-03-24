# purpose - make Figs. 3a and 3b showing CoPRO read2 (TSS) and read1 (paused PolII) and BNase-seq, all three with all insert-sized reads.

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

#set bedfile of plus 1 nucleosome, sorted by decreasing RNA-expression
Plus1_BEDFILE_expressed=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed
Plus1_BEDFILE_unexpressed=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230720_plus1_minus1/output_v2_NonRed_OCT_Hex_Tet_Unexpressed_230905/K562_Plus1_NOsort_Unexpressed_nonRedOct_Hex_Tet.bed

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig3_atTSS/fig3a_b_v2_230926_output

#set bam library file to BI_rep1
BAM1=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/MEP_Project/01_BAM/CoPRO/CoPRO_K562_MERGE.bam
BAM2=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/K562_benzonase-seq_master.bam
BAM3=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_ChIPs/MERGED_datasets/19354_19355_NoBenz_10sonicCycles_XO_PolII_master.bam

#set blacklist
BLACKLIST=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/hg19_Blacklist.bed

#set scriptmanager and job
SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

mkdir -p $OUTPUT

JOBSTATS="#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=24gb
#PBS -l walltime=10:00:00
#PBS -A open

source ~/.bashrc #configures shell to use conda activate
conda activate bioinfo"

#set output file names
Plus1_BEDFILE=$(echo $Plus1_BEDFILE_expressed | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_expressedANDunexpressed.bed"}')
Plus1_BEDFILE_2000=$(echo $Plus1_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2000bp.bed"}')
Plus1_BEDFILE_3000=$(echo $Plus1_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_3000bp.bed"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2.out"}')
CDT1=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2"}')
CDT1_sense_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2_sense.cdt.gz"}')
CDT1_anti_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2_anti.cdt.gz"}')
CDT1_sense=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2_sense.cdt"}')
CDT1_anti=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2_anti.cdt"}')
CDT1_SCALED_sense=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2_sense_scaled.cdt"}')
CDT1_SCALED_anti=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read2_anti_scaled.cdt"}')
SCALE1=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_read2"}')
SCALE1a=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_read2_ScalingFactors.out"}')
PNG1_sense=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read2_sense.png"}')
PNG1_anti=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read2_anti.png"}')
PNG1_merge=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read2_merged.png"}')
SVG1=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read2_labeled.svg"}')
OUT1a_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read2.out"}')
CDT1d_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read2"}')
CDT1_sense_3000_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read2_sense.cdt.gz"}')
CDT1_anti_3000_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read2_anti.cdt.gz"}')
CDT1_sense_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read2_sense.cdt"}')
CDT1_anti_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read2_anti.cdt"}')
SCALE1_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite_read2"}')
SCALE1a_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite_read2_ScalingFactors.out"}')
CDT1_SCALED_COMP_sense_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read2_sense.cdt"}')
CDT1_SCALED_COMP_anti_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read2_anti.cdt"}')
SCALED_OUT1_sense_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read2_sense.tab"}')
SCALED_OUT1_anti_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read2_anti.tab"}')
OUT2=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1.out"}')
CDT2=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1"}')
CDT2_sense_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1_sense.cdt.gz"}')
CDT2_anti_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1_anti.cdt.gz"}')
CDT2_sense=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1_sense.cdt"}')
CDT2_anti=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1_anti.cdt"}')
CDT2_SCALED_sense=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1_sense_scaled.cdt"}')
CDT2_SCALED_anti=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_read1_anti_scaled.cdt"}')
SCALE2=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_read1"}')
SCALE2a=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_read1_ScalingFactors.out"}')
PNG2_sense=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read1_sense.png"}')
PNG2_anti=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read1_anti.png"}')
PNG2_merge=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read1_merged.png"}')
SVG2=$(echo "$BAM1a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_read1_labeled.svg"}')
OUT2a_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read1.out"}')
CDT2d_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read1"}')
CDT2_sense_3000_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read1_sense.cdt.gz"}')
CDT2_anti_3000_gz=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read1_anti.cdt.gz"}')
CDT2_sense_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read1_sense.cdt"}')
CDT2_anti_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_read1_anti.cdt"}')
SCALE2_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite_read1"}')
SCALE2a_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite_read1_ScalingFactors.out"}')
CDT2_SCALED_COMP_sense_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.cdt"}')
CDT2_SCALED_COMP_anti_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.cdt"}')
SCALED_OUT2_sense_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read1_sense.tab"}')
SCALED_OUT2_anti_3000=$(echo "$BAM1a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_read1_anti.tab"}')
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT3=$(echo "$BAM2a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT3=$(echo "$BAM2a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1""}')
CDT3b=$(echo "$BAM2a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT3c=$(echo "$BAM2a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT3_SCALED=$(echo "$BAM2a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE3=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT"}')
SCALE3a=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG3=$(echo "$BAM2a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG3=$(echo "$BAM2a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT3a_3000=$(echo "$BAM2a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT3d_3000=$(echo "$BAM2a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT3e_3000=$(echo "$BAM2a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT3f_3000=$(echo "$BAM2a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
SCALE3_OUT=$(echo "$BAM2a" | awk -F. '{print $1"_ForComposite"}')
SCALE3a_OUT=$(echo "$BAM2a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT3_SCALED_COMP_3000=$(echo "$BAM2a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT3_3000=$(echo "$BAM2a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT4=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"read1.out"}')
CDT4=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1}')
CDT4_sense_gz=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_sense.cdt.gz"}')
CDT4_anti_gz=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_anti.cdt.gz"}')
CDT4_sense=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_sense.cdt"}')
CDT4_anti=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_anti.cdt"}')
CDT4_SCALED_sense=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_sense_scaled.cdt"}')
CDT4_SCALED_anti=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_anti_scaled.cdt"}')
SCALE4=$(echo "$BAM3a" | awk -F. '{print $1"_ForCDT"}')
SCALE4a=$(echo "$BAM3a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG4_sense=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_sense.png"}')
PNG4_anti=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_anti.png"}')
PNG4_merge=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_merged.png"}')
SVG4=$(echo "$BAM3a""_""$Plus1_BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT4a_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite.out"}')
CDT4d_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT4_sense_3000_gz=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_sense.cdt.gz"}')
CDT4_anti_3000_gz=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_anti.cdt.gz"}')
CDT4_sense_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_sense.cdt"}')
CDT4_anti_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_anti.cdt"}')
SCALE4_OUT=$(echo "$BAM3a" | awk -F. '{print $1"_ForComposite"}')
SCALE4a_OUT=$(echo "$BAM3a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT4_SCALED_COMP_sense_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_sense.cdt"}')
CDT4_SCALED_COMP_anti_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_anti.cdt"}')
SCALED_OUT4_sense_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_sense.tab"}')
SCALED_OUT4_anti_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_anti.tab"}')
SCALED_OUT4_final_3000=$(echo "$BAM3a""_""$Plus1_BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled_final.tab"}')

sampleID=fig3a_b_v2_separateStrands.pbs
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#concatenate TSS bedfiles (expressed then unexpressed)" >> $sampleID
echo "cat $Plus1_BEDFILE_expressed $Plus1_BEDFILE_unexpressed > $Plus1_BEDFILE" >> $sampleID
echo "#expand bedfiles by 2000 bp and 3000 bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2000 $Plus1_BEDFILE -o=$Plus1_BEDFILE_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=3000 $Plus1_BEDFILE -o=$Plus1_BEDFILE_3000" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 -z --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $Plus1_BEDFILE_2000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $Plus1_BEDFILE_2000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $Plus1_BEDFILE_2000 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT4 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 $Plus1_BEDFILE_2000 $BAM3" >> $sampleID
echo "#scale output files: options - total tag scaling -t" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE4 $BAM3" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1_sense_gz > $CDT1_sense" >> $sampleID
echo "gunzip -c $CDT1_anti_gz > $CDT1_anti" >> $sampleID
echo "gunzip -c $CDT2_sense_gz > $CDT2_sense" >> $sampleID
echo "gunzip -c $CDT2_anti_gz > $CDT2_anti" >> $sampleID
echo "gunzip -c $CDT3b > $OUTPUT/$CDT3c" >> $sampleID
echo "gunzip -c $CDT4_sense_gz > $CDT4_sense" >> $sampleID
echo "gunzip -c $CDT4_anti_gz > $CDT4_anti" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_sense --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_anti --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_sense --scaling-factor=\$(cat $SCALE2a | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_anti --scaling-factor=\$(cat $SCALE2a | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT3c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_sense --scaling-factor=\$(cat $SCALE4a | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_anti --scaling-factor=\$(cat $SCALE4a | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_anti" >> $sampleID
echo "#make heatmap **change threshold if necessary" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG1_sense -c ff0000 -p=0.95 $CDT1_SCALED_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG1_anti -c 0000ff -p=0.95 $CDT1_SCALED_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG2_sense -c 7f00ff -p=0.95 $CDT2_SCALED_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG2_anti -c 0099ff -p=0.95 $CDT2_SCALED_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG3 -c ff00ff -p=0.95 $CDT3_SCALED" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG4_sense -c 833c0c -p=0.95 $CDT4_SCALED_sense" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG4_anti -c c65911 -p=0.95 $CDT4_SCALED_anti" >> $sampleID
echo "#merge heatmaps" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation merge-heatmap -o=$PNG1_merge $PNG1_sense $PNG1_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation merge-heatmap -o=$PNG2_merge $PNG2_sense $PNG2_anti" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation merge-heatmap -o=$PNG4_merge $PNG4_sense $PNG4_anti" >> $sampleID
echo "#label above heatmaps" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG1_merge --output=$SVG1 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='18,236 Nucleosome positions'" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG2_merge --output=$SVG2 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='18,236 Nucleosome positions'" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG3 --output=$SVG3 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='18,236 Nucleosome positions'" >> $sampleID
echo "java -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG4_merge --output=$SVG4 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='18,236 Nucleosome positions'" >> $sampleID
echo "#prep files for plotter" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -2 -z --output-matrix=$CDT1d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1a_3000 $Plus1_BEDFILE_3000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT2d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2a_3000 $Plus1_BEDFILE_3000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3a_3000 $Plus1_BEDFILE_3000 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --output-matrix=$CDT4d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4a_3000 $Plus1_BEDFILE_3000 $BAM3" >> $sampleID
echo "#scale output files: options - total tag scaling -t; *this is same regardless of bedfile" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1_OUT $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2_OUT $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3_OUT $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE4_OUT $BAM3" >> $sampleID
echo "#unzip cdt files" >> $sampleID
echo "gunzip -c $CDT1_sense_3000_gz > $CDT1_sense_3000" >> $sampleID
echo "gunzip -c $CDT1_anti_3000_gz > $CDT1_anti_3000" >> $sampleID
echo "gunzip -c $CDT2_sense_3000_gz > $CDT2_sense_3000" >> $sampleID
echo "gunzip -c $CDT2_anti_3000_gz > $CDT2_anti_3000" >> $sampleID
echo "gunzip -c $CDT3e_3000 > $OUTPUT/$CDT3f_3000" >> $sampleID
echo "gunzip -c $CDT4_sense_3000_gz > $CDT4_sense_3000" >> $sampleID
echo "gunzip -c $CDT4_anti_3000_gz > $CDT4_anti_3000" >> $sampleID
echo "#scale CDT files data in matrix by scaling factor" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_COMP_sense_3000 --scaling-factor=\$(cat $SCALE1a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_sense_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_COMP_anti_3000 --scaling-factor=\$(cat $SCALE1a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT1_anti_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_COMP_sense_3000 --scaling-factor=\$(cat $SCALE2a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_sense_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_COMP_anti_3000 --scaling-factor=\$(cat $SCALE2a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT2_anti_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE3a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT3f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_sense_3000 --scaling-factor=\$(cat $SCALE4a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_sense_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_anti_3000 --scaling-factor=\$(cat $SCALE4a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT4_anti_3000" >> $sampleID
echo "#make scaled OUT file for each strand" >> $sampleID
echo "perl $JOB $CDT1_SCALED_COMP_sense_3000 $SCALED_OUT1_sense_3000" >> $sampleID
echo "perl $JOB $CDT1_SCALED_COMP_anti_3000 $SCALED_OUT1_anti_3000" >> $sampleID
echo "perl $JOB $CDT2_SCALED_COMP_sense_3000 $SCALED_OUT2_sense_3000" >> $sampleID
echo "perl $JOB $CDT2_SCALED_COMP_anti_3000 $SCALED_OUT2_anti_3000" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_3000 $SCALED_OUT3_3000" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_sense_3000 $SCALED_OUT4_sense_3000" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_anti_3000 $SCALED_OUT4_anti_3000" >> $sampleID
echo "#concatenate OUT fles and take lines 1,2,4 to final composite files for each library. #Final version not made for subfigures 3 as above files "SCALED_OUT_3000" are their final form." >> $sampleID
echo "cat $SCALED_OUT4_sense_3000 $SCALED_OUT4_anti_3000 | awk 'NR==1;NR==2;NR==4' > $SCALED_OUT4_final_3000" >> $sampleID
echo "#remove intermediate files"
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1_sense_gz" >> $sampleID
echo "rm $CDT1_anti_gz" >> $sampleID
echo "rm $CDT1_sense" >> $sampleID
echo "rm $CDT1_anti" >> $sampleID
echo "rm $CDT1_SCALED_sense" >> $sampleID
echo "rm $CDT1_SCALED_anti" >> $sampleID
echo "rm $SCALE1a" >> $sampleID
echo "rm $PNG1_sense" >> $sampleID
echo "rm $PNG1_anti" >> $sampleID
echo "rm $PNG1_merge" >> $sampleID
echo "rm $OUT1a_3000" >> $sampleID
echo "rm $CDT1_sense_3000_gz" >> $sampleID
echo "rm $CDT1_anti_3000_gz" >> $sampleID
echo "rm $CDT1_sense_3000" >> $sampleID
echo "rm $CDT1_anti_3000" >> $sampleID
echo "rm $SCALE1a_OUT" >> $sampleID
echo "rm $CDT1_SCALED_COMP_sense_3000" >> $sampleID
echo "rm $CDT1_SCALED_COMP_anti_3000" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2_sense_gz" >> $sampleID
echo "rm $CDT2_anti_gz" >> $sampleID
echo "rm $CDT2_sense" >> $sampleID
echo "rm $CDT2_anti" >> $sampleID
echo "rm $CDT2_SCALED_sense" >> $sampleID
echo "rm $CDT2_SCALED_anti" >> $sampleID
echo "rm $SCALE2a" >> $sampleID
echo "rm $PNG2_sense" >> $sampleID
echo "rm $PNG2_anti" >> $sampleID
echo "rm $PNG2_merge" >> $sampleID
echo "rm $OUT2a_3000" >> $sampleID
echo "rm $CDT2_sense_3000_gz" >> $sampleID
echo "rm $CDT2_anti_3000_gz" >> $sampleID
echo "rm $CDT2_sense_3000" >> $sampleID
echo "rm $CDT2_anti_3000" >> $sampleID
echo "rm $SCALE2a_OUT" >> $sampleID
echo "rm $CDT2_SCALED_COMP_sense_3000" >> $sampleID
echo "rm $CDT2_SCALED_COMP_anti_3000" >> $sampleID
echo "rm $OUT3" >> $sampleID
echo "rm $CDT3b" >> $sampleID
echo "rm $CDT3c" >> $sampleID
echo "rm $CDT3_SCALED" >> $sampleID
echo "rm $SCALE3a" >> $sampleID
echo "rm $PNG3" >> $sampleID
echo "rm $OUT3a_3000" >> $sampleID
echo "rm $CDT3e_3000" >> $sampleID
echo "rm $CDT3f_3000" >> $sampleID
echo "rm $SCALE3a_OUT" >> $sampleID
echo "rm $CDT3_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT4" >> $sampleID
echo "rm $CDT4_sense_gz" >> $sampleID
echo "rm $CDT4_anti_gz" >> $sampleID
echo "rm $CDT4_sense" >> $sampleID
echo "rm $CDT4_anti" >> $sampleID
echo "rm $CDT4_SCALED_sense" >> $sampleID
echo "rm $CDT4_SCALED_anti" >> $sampleID
echo "rm $SCALE4a" >> $sampleID
echo "rm $PNG4_sense" >> $sampleID
echo "rm $PNG4_anti" >> $sampleID
echo "rm $PNG4_merge" >> $sampleID
echo "rm $OUT4a_3000" >> $sampleID
echo "rm $CDT4_sense_3000_gz" >> $sampleID
echo "rm $CDT4_anti_3000_gz" >> $sampleID
echo "rm $CDT4_sense_3000" >> $sampleID
echo "rm $CDT4_anti_3000" >> $sampleID
echo "rm $SCALE4a_OUT" >> $sampleID
echo "rm $CDT4_SCALED_COMP_sense_3000" >> $sampleID
echo "rm $CDT4_SCALED_COMP_anti_3000" >> $sampleID
echo "rm $SCALED_OUT4_sense_3000" >> $sampleID
echo "rm $SCALED_OUT4_anti_3000" >> $sampleID
echo "# finish script" >> $sampleID
