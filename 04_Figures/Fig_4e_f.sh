# purpose - make Fig. 3e/f showing Cut & Run of various targets at the +1 Nucleosome dyad with size selection of 128-164 bp. This version (v3) uses all sites for composite plots and uses size selection for heatmaps and tag pile-ups.

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

#set bedfiles with Plus1 (downstream) nucleosome positions based on non-redundant bedfile-based calls
BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230720_plus1_minus1/output_v2_NonRed_Oct_Hex_Tet_230825/K562_Plus1_SORTbyRNAexp_nonRedOct_Hex_Tet.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/nre_fig3_240219/fig3e_f_CUTandRUN_128to164_output_v3_240318

#set bam library file to BI_rep1
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN/H2AZ/H2AZ_4DN_CUTRUN.bam
BAM2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN/H3K4me1/H3K4me1_4DN_CUTRUN.bam
BAM3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN/H3K4me3/H3K4me3_4DN_CUTRUN.bam
BAM4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN/H3K27ac/H3K27ac_4DN_CUTRUN.bam
BAM5=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN/H3K27me3/H3K27me3_4DN_CUTRUN.bam
BAM6=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN/CTCF/CTCF_4DN_CUTRUN.bam
BAM7=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/00_BAM/CUTRUN/IgG/IgG_4DN_CUTRUN.bam

#set blacklist
BLACKLIST=/storage/group/bfp2/default/juk398-JordanKrebs/hg19_Blacklist.bed

#set scriptmanager and job
SCRIPTMANAGER=/storage/group/bfp2/default/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/jobs/sum_Col_CDT.pl

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
#SBATCH --time=6:00:00
#SBATCH --partition=open

source ~/.bashrc #configures shell to use conda activate
module load anaconda
conda activate bioinfo"

#set output file names
BEDFILE_2000=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_2000bp.bed"}')
BEDFILE_3000=$(echo $BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1"_3000bp.bed"}')
BAM1a=$(echo "$BAM1" | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$BEDFILE_2000" | awk -F. '{print $1}')
CDT1b=$(echo "$BAM1a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT1c=$(echo "$BAM1a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT1_SCALED=$(echo "$BAM1a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE1=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT"}')
SCALE1a=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG1=$(echo "$BAM1a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG1=$(echo "$BAM1a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT2=$(echo "$BAM2a""_""$BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT2=$(echo "$BAM2a""_""$BEDFILE_2000" | awk -F. '{print $1}')
CDT2b=$(echo "$BAM2a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT2c=$(echo "$BAM2a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT2_SCALED=$(echo "$BAM2a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE2=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT"}')
SCALE2a=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG2=$(echo "$BAM2a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG2=$(echo "$BAM2a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT3=$(echo "$BAM3a""_""$BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT3=$(echo "$BAM3a""_""$BEDFILE_2000" | awk -F. '{print $1}')
CDT3b=$(echo "$BAM3a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT3c=$(echo "$BAM3a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT3_SCALED=$(echo "$BAM3a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE3=$(echo "$BAM3a" | awk -F. '{print $1"_ForCDT"}')
SCALE3a=$(echo "$BAM3a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG3=$(echo "$BAM3a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG3=$(echo "$BAM3a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
BAM4a=$(echo $BAM4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT4=$(echo "$BAM4a""_""$BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT4=$(echo "$BAM4a""_""$BEDFILE_2000" | awk -F. '{print $1}')
CDT4b=$(echo "$BAM4a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT4c=$(echo "$BAM4a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT4_SCALED=$(echo "$BAM4a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE4=$(echo "$BAM4a" | awk -F. '{print $1"_ForCDT"}')
SCALE4a=$(echo "$BAM4a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG4=$(echo "$BAM4a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG4=$(echo "$BAM4a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
BAM5a=$(echo $BAM5 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT5=$(echo "$BAM5a""_""$BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT5=$(echo "$BAM5a""_""$BEDFILE_2000" | awk -F. '{print $1}')
CDT5b=$(echo "$BAM5a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT5c=$(echo "$BAM5a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT5_SCALED=$(echo "$BAM5a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE5=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT"}')
SCALE5a=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG5=$(echo "$BAM5a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG5=$(echo "$BAM5a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
BAM6a=$(echo $BAM6 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT6=$(echo "$BAM6a""_""$BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT6=$(echo "$BAM6a""_""$BEDFILE_2000" | awk -F. '{print $1}')
CDT6b=$(echo "$BAM6a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT6c=$(echo "$BAM6a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT6_SCALED=$(echo "$BAM6a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE6=$(echo "$BAM6a" | awk -F. '{print $1"_ForCDT"}')
SCALE6a=$(echo "$BAM6a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG6=$(echo "$BAM6a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG6=$(echo "$BAM6a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
BAM7a=$(echo $BAM7 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT7=$(echo "$BAM7a""_""$BEDFILE_2000" | awk -F. '{print $1"_midpoint.out"}')
CDT7=$(echo "$BAM7a""_""$BEDFILE_2000" | awk -F. '{print $1}')
CDT7b=$(echo "$BAM7a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT7c=$(echo "$BAM7a""_""$BEDFILE_2000" | awk -F. '{print $1"_combined.cdt"}')
CDT7_SCALED=$(echo "$BAM7a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.cdt"}')
SCALE7=$(echo "$BAM7a" | awk -F. '{print $1"_ForCDT"}')
SCALE7a=$(echo "$BAM7a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG7=$(echo "$BAM7a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled.png"}')
SVG7=$(echo "$BAM7a""_""$BEDFILE_2000" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT1a_3000=$(echo "$BAM1a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT1d_3000=$(echo "$BAM1a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT1e_3000=$(echo "$BAM1a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT1f_3000=$(echo "$BAM1a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT1_SCALED_COMP_3000=$(echo "$BAM1a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT1_3000=$(echo "$BAM1a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT2a_3000=$(echo "$BAM2a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT2d_3000=$(echo "$BAM2a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT2e_3000=$(echo "$BAM2a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT2f_3000=$(echo "$BAM2a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT2_SCALED_COMP_3000=$(echo "$BAM2a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT2_3000=$(echo "$BAM2a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT3a_3000=$(echo "$BAM3a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT3d_3000=$(echo "$BAM3a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT3e_3000=$(echo "$BAM3a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT3f_3000=$(echo "$BAM3a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT3_SCALED_COMP_3000=$(echo "$BAM3a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT3_3000=$(echo "$BAM3a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT4a_3000=$(echo "$BAM4a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT4d_3000=$(echo "$BAM4a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT4e_3000=$(echo "$BAM4a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT4f_3000=$(echo "$BAM4a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT4_SCALED_COMP_3000=$(echo "$BAM4a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT4_3000=$(echo "$BAM4a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT5a_3000=$(echo "$BAM5a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT5d_3000=$(echo "$BAM5a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT5e_3000=$(echo "$BAM5a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT5f_3000=$(echo "$BAM5a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT5_SCALED_COMP_3000=$(echo "$BAM5a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT5_3000=$(echo "$BAM5a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT6a_3000=$(echo "$BAM6a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT6d_3000=$(echo "$BAM6a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT6e_3000=$(echo "$BAM6a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT6f_3000=$(echo "$BAM6a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT6_SCALED_COMP_3000=$(echo "$BAM6a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT6_3000=$(echo "$BAM6a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
OUT7a_3000=$(echo "$BAM7a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT7d_3000=$(echo "$BAM7a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite"}')
CDT7e_3000=$(echo "$BAM7a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT7f_3000=$(echo "$BAM7a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT7_SCALED_COMP_3000=$(echo "$BAM7a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT7_3000=$(echo "$BAM7a""_""$BEDFILE_3000" | awk -F. '{print $1"_ForComposite_scaled.tab"}')

sampleID=fig3e_f_CUTandRUN_+1_128to164bp_v3_240318.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#expand bedfiles by 2000 bp and 3000 bp" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2000 $BEDFILE -o=$BEDFILE_2000" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=3000 $BEDFILE -o=$BEDFILE_3000" >> $sampleID
echo "#scale output files: options - total tag scaling -t" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE4 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE5 $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE6 $BAM6" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE7 $BAM7" >> $sampleID
echo "#do analysis for Nucleosome position" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 --min-insert=128 --max-insert=164 $BEDFILE_2000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 --min-insert=128 --max-insert=164 $BEDFILE_2000 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 --min-insert=128 --max-insert=164 $BEDFILE_2000 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 --min-insert=128 --max-insert=164 $BEDFILE_2000 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 --min-insert=128 --max-insert=164 $BEDFILE_2000 $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT6 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6 --min-insert=128 --max-insert=164 $BEDFILE_2000 $BAM6" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT7 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT7 --min-insert=128 --max-insert=164 $BEDFILE_2000 $BAM7" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "gunzip -c $CDT1b > $CDT1c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1c" >> $sampleID
echo "gunzip -c $CDT2b > $CDT2c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED --scaling-factor=\$(cat $SCALE2a | cut -f2 | tail -1 | awk '{print \$1}') $CDT2c" >> $sampleID
echo "gunzip -c $CDT3b > $CDT3c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT3c" >> $sampleID
echo "gunzip -c $CDT4b > $CDT4c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED --scaling-factor=\$(cat $SCALE4a | cut -f2 | tail -1 | awk '{print \$1}') $CDT4c" >> $sampleID
echo "gunzip -c $CDT5b > $CDT5c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT5c" >> $sampleID
echo "gunzip -c $CDT6b > $CDT6c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6c" >> $sampleID
echo "gunzip -c $CDT7b > $CDT7c" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT7_SCALED --scaling-factor=\$(cat $SCALE7a | cut -f2 | tail -1 | awk '{print \$1}') $CDT7c" >> $sampleID
echo "#make heatmaps, merge heatmaps [colors set] **change threshold if necessary" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG1 -c 999900 -p=0.95 $CDT1_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG2 -c ff6699 -p=0.95 $CDT2_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG3 -c ff9900 -p=0.95 $CDT3_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG4 -c 9900ff -p=0.95 $CDT4_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG5 -c 999999 -p=0.95 $CDT5_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG6 -c 703000 -p=0.95 $CDT6_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG7 -c 000000 -p=0.95 $CDT7_SCALED" >> $sampleID
echo "#label above heatmaps" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG1 --output=$SVG1 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='11,714 Nucleosome positions'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG2 --output=$SVG2 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='11,714 Nucleosome positions'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG3 --output=$SVG3 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='11,714 Nucleosome positions'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG4 --output=$SVG4 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='11,714 Nucleosome positions'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG5 --output=$SVG5 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='11,714 Nucleosome positions'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG6 --output=$SVG6 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='11,714 Nucleosome positions'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG7 --output=$SVG7 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from +1 Dyad (kb)' --y-label='11,714 Nucleosome positions'" >> $sampleID
echo "#prep files for plotter" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1a_3000 --min-insert=128 --max-insert=164 $BEDFILE_3000 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2a_3000 --min-insert=128 --max-insert=164 $BEDFILE_3000 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3a_3000 --min-insert=128 --max-insert=164 $BEDFILE_3000 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4a_3000 --min-insert=128 --max-insert=164 $BEDFILE_3000 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5a_3000 --min-insert=128 --max-insert=164 $BEDFILE_3000 $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT6d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT6a_3000 --min-insert=128 --max-insert=164 $BEDFILE_3000 $BAM6" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT7d_3000 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT7a_3000 --min-insert=128 --max-insert=164 $BEDFILE_3000 $BAM7" >> $sampleID
echo "#scale data in matrix by scaling factor" >> $sampleID
echo "gunzip -c $CDT1e_3000 > $OUTPUT/$CDT1f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE1a | cut -f2 | tail -1 | awk '{print \$1}') $CDT1f_3000" >> $sampleID
echo "gunzip -c $CDT2e_3000 > $OUTPUT/$CDT2f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE2a | cut -f2 | tail -1 | awk '{print \$1}') $CDT2f_3000" >> $sampleID
echo "gunzip -c $CDT3e_3000 > $OUTPUT/$CDT3f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE3a | cut -f2 | tail -1 | awk '{print \$1}') $CDT3f_3000" >> $sampleID
echo "gunzip -c $CDT4e_3000 > $OUTPUT/$CDT4f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE4a | cut -f2 | tail -1 | awk '{print \$1}') $CDT4f_3000" >> $sampleID
echo "gunzip -c $CDT5e_3000 > $OUTPUT/$CDT5f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE5a | cut -f2 | tail -1 | awk '{print \$1}') $CDT5f_3000" >> $sampleID
echo "gunzip -c $CDT6e_3000 > $OUTPUT/$CDT6f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT6_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE6a | cut -f2 | tail -1 | awk '{print \$1}') $CDT6f_3000" >> $sampleID
echo "gunzip -c $CDT7e_3000 > $OUTPUT/$CDT7f_3000" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT7_SCALED_COMP_3000 --scaling-factor=\$(cat $SCALE7a | cut -f2 | tail -1 | awk '{print \$1}') $CDT7f_3000" >> $sampleID
echo "#make scale OUT file" >> $sampleID
echo "perl $JOB $CDT1_SCALED_COMP_3000 $SCALED_OUT1_3000" >> $sampleID
echo "perl $JOB $CDT2_SCALED_COMP_3000 $SCALED_OUT2_3000" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_3000 $SCALED_OUT3_3000" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_3000 $SCALED_OUT4_3000" >> $sampleID
echo "perl $JOB $CDT5_SCALED_COMP_3000 $SCALED_OUT5_3000" >> $sampleID
echo "perl $JOB $CDT6_SCALED_COMP_3000 $SCALED_OUT6_3000" >> $sampleID
echo "perl $JOB $CDT7_SCALED_COMP_3000 $SCALED_OUT7_3000" >> $sampleID
echo "#remove intermediate files"
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1b" >> $sampleID
echo "rm $CDT1c" >> $sampleID
echo "rm $CDT1_SCALED" >> $sampleID
echo "rm $SCALE1a" >> $sampleID
echo "rm $PNG1" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2b" >> $sampleID
echo "rm $CDT2c" >> $sampleID
echo "rm $CDT2_SCALED" >> $sampleID
echo "rm $SCALE2a" >> $sampleID
echo "rm $PNG2" >> $sampleID
echo "rm $OUT3" >> $sampleID
echo "rm $CDT3b" >> $sampleID
echo "rm $CDT3c" >> $sampleID
echo "rm $CDT3_SCALED" >> $sampleID
echo "rm $SCALE3a" >> $sampleID
echo "rm $PNG3" >> $sampleID
echo "rm $OUT4" >> $sampleID
echo "rm $CDT4b" >> $sampleID
echo "rm $CDT4c" >> $sampleID
echo "rm $CDT4_SCALED" >> $sampleID
echo "rm $SCALE4a" >> $sampleID
echo "rm $PNG4" >> $sampleID
echo "rm $OUT5" >> $sampleID
echo "rm $CDT5b" >> $sampleID
echo "rm $CDT5c" >> $sampleID
echo "rm $CDT5_SCALED" >> $sampleID
echo "rm $SCALE5a" >> $sampleID
echo "rm $PNG5" >> $sampleID
echo "rm $OUT6" >> $sampleID
echo "rm $CDT6b" >> $sampleID
echo "rm $CDT6c" >> $sampleID
echo "rm $CDT6_SCALED" >> $sampleID
echo "rm $SCALE6a" >> $sampleID
echo "rm $PNG6" >> $sampleID
echo "rm $OUT7" >> $sampleID
echo "rm $CDT7b" >> $sampleID
echo "rm $CDT7c" >> $sampleID
echo "rm $CDT7_SCALED" >> $sampleID
echo "rm $SCALE7a" >> $sampleID
echo "rm $PNG7" >> $sampleID
echo "rm $OUT1a_3000" >> $sampleID
echo "rm $CDT1e_3000" >> $sampleID
echo "rm $CDT1f_3000" >> $sampleID
echo "rm $CDT1_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT2a_3000" >> $sampleID
echo "rm $CDT2e_3000" >> $sampleID
echo "rm $CDT2f_3000" >> $sampleID
echo "rm $CDT2_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT3a_3000" >> $sampleID
echo "rm $CDT3e_3000" >> $sampleID
echo "rm $CDT3f_3000" >> $sampleID
echo "rm $CDT3_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT4a_3000" >> $sampleID
echo "rm $CDT4e_3000" >> $sampleID
echo "rm $CDT4f_3000" >> $sampleID
echo "rm $CDT4_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT5a_3000" >> $sampleID
echo "rm $CDT5e_3000" >> $sampleID
echo "rm $CDT5f_3000" >> $sampleID
echo "rm $CDT5_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT6a_3000" >> $sampleID
echo "rm $CDT6e_3000" >> $sampleID
echo "rm $CDT6f_3000" >> $sampleID
echo "rm $CDT6_SCALED_COMP_3000" >> $sampleID
echo "rm $OUT7a_3000" >> $sampleID
echo "rm $CDT7e_3000" >> $sampleID
echo "rm $CDT7f_3000" >> $sampleID
echo "rm $CDT7_SCALED_COMP_3000" >> $sampleID
echo "# finish script" >> $sampleID
