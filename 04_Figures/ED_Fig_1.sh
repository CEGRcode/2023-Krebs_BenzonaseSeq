# purpose - make supplemtnal Fig. 1 showing benzonase-seq and MNase-seq at ref-seq TSS(s) sorted by K562 gene expression.

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
TSS_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/bedfiles/K562_CoPRO-expressed_Gene-refSeqTSS_2000bp.bed
CpG_BEDFILE=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/bedfiles/UCSCgb_hg19_CpGislands_230426.bed

#output directory
OUTPUT=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/figures/fig1_atTSS_CpGsort/ExtData1_240102_output

#set bam library file to BI_rep1
BAM1=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230718_MERGE/K562_benzonase-seq_master.bam
BAM2=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_MNase_DNase/final_files/SRR3211679_master.bam
BAM3=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_MNase_DNase/final_files/SRR3211680_master.bam
BAM4=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_MNase_DNase/final_files/SRR3211681_master.bam
BAM5=/storage/group/bfp2/default/juk398-JordanKrebs/NucleosomeAtlas_project/230810_MNase_DNase/final_files/SRR3211682_master.bam

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
TSS_BEDFILEa=$(echo $TSS_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
CpG_BEDFILEa=$(echo $CpG_BEDFILE | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
BEDFILE=K562_CoPRO_TSS_sortByCpGIsland.bed
BEDFILEb=K562_CoPRO_TSS_sortByCpGIsland_2200bp.bed
ALIGN_BED1=$(echo "$CpG_BEDFILEa""_""$TSS_BEDFILEa" | awk -F. '{print $1".cdt"}')
ALIGN_BED2=$(echo "$CpG_BEDFILEa""_K562_CoPRO_TSS_sortByCpGIsland" | awk -F. '{print $1".cdt"}')
CpG_PNG=$(echo "$CpG_BEDFILEa""_K562_CoPRO_TSS_sortByCpGIsland" | awk -F. '{print $1".png"}')
CpG_SVG=$(echo "$CpG_BEDFILEa""_K562_CoPRO_TSS_sortByCpGIsland" | awk -F. '{print $1".svg"}')
TOP_BED=K562_CoPRO_TSS_sortByCpGIsland_top2500.bed
BOTTOM_BED=K562_CoPRO_TSS_sortByCpGIsland_bottom1600.bed
TOP_BED_2200=K562_CoPRO_TSS_sortByCpGIsland_top2500_2200bp.bed
BOTTOM_BED_2200=K562_CoPRO_TSS_sortByCpGIsland_bottom1600_2200.bed
ALIGN_BED_top=$(echo "K562_CoPRO_TSS_sortByCpGIsland_top2500_2200bp" | awk -F. '{print $1".cdt"}')
ALIGN_BED_bottom=$(echo "K562_CoPRO_TSS_sortByCpGIsland_bottom1600_2200bp" | awk -F. '{print $1".cdt"}')
OUT_top=$(echo "K562_CoPRO_TSS_sortByCpGIsland_top2500_2200bp" | awk -F. '{print $1".out"}')
OUT_bottom=$(echo "K562_CoPRO_TSS_sortByCpGIsland_bottom1600_2200bp" | awk -F. '{print $1".out"}')
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT1=$(echo "$BAM1a""_""$BEDFILE" | awk -F. '{print $1"_midpoint.out"}')
CDT1=$(echo "$BAM1a""_""$BEDFILE" | awk -F. '{print $1}')
CDT1b=$(echo "$BAM1a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT1c=$(echo "$BAM1a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt"}')
CDT1_SCALED=$(echo "$BAM1a""_""$BEDFILE" | awk -F. '{print $1"_scaled.cdt"}')
SCALE1=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT"}')
SCALE1a=$(echo "$BAM1a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG1=$(echo "$BAM1a""_""$BEDFILE" | awk -F. '{print $1"_scaled.png"}')
SVG1=$(echo "$BAM1a""_""$BEDFILE" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT1a_top=$(echo "$BAM1a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
OUT1a_bottom=$(echo "$BAM1a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT1d_top=$(echo "$BAM1a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT1d_bottom=$(echo "$BAM1a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT1e_top=$(echo "$BAM1a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT1e_bottom=$(echo "$BAM1a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT1f_top=$(echo "$BAM1a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT1f_bottom=$(echo "$BAM1a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
SCALE1_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite"}')
SCALE1a_OUT=$(echo "$BAM1a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT1_SCALED_COMP_top=$(echo "$BAM1a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
CDT1_SCALED_COMP_bottom=$(echo "$BAM1a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT1_top=$(echo "$BAM1a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
SCALED_OUT1_bottom=$(echo "$BAM1a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_midpoint.out"}')
CDT2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1}')
CDT2b=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT2c=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt"}')
CDT2_SCALED=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_scaled.cdt"}')
SCALE2=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT"}')
SCALE2a=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_scaled.png"}')
SVG2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT2a_top=$(echo "$BAM2a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
OUT2a_bottom=$(echo "$BAM2a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT2d_top=$(echo "$BAM2a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT2d_bottom=$(echo "$BAM2a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT2e_top=$(echo "$BAM2a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT2e_bottom=$(echo "$BAM2a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT2f_top=$(echo "$BAM2a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT2f_bottom=$(echo "$BAM2a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
SCALE2_OUT=$(echo "$BAM2a" | awk -F. '{print $1"_ForComposite"}')
SCALE2a_OUT=$(echo "$BAM2a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT2_SCALED_COMP_top=$(echo "$BAM2a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
CDT2_SCALED_COMP_bottom=$(echo "$BAM2a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT2_top=$(echo "$BAM2a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
SCALED_OUT2_bottom=$(echo "$BAM2a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT3=$(echo "$BAM3a""_""$BEDFILE" | awk -F. '{print $1"_midpoint.out"}')
CDT3=$(echo "$BAM3a""_""$BEDFILE" | awk -F. '{print $1}')
CDT3b=$(echo "$BAM3a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT3c=$(echo "$BAM3a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt"}')
CDT3_SCALED=$(echo "$BAM3a""_""$BEDFILE" | awk -F. '{print $1"_scaled.cdt"}')
SCALE3=$(echo "$BAM3a" | awk -F. '{print $1"_ForCDT"}')
SCALE3a=$(echo "$BAM3a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG3=$(echo "$BAM3a""_""$BEDFILE" | awk -F. '{print $1"_scaled.png"}')
SVG3=$(echo "$BAM3a""_""$BEDFILE" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT3a_top=$(echo "$BAM3a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
OUT3a_bottom=$(echo "$BAM3a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT3d_top=$(echo "$BAM3a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT3d_bottom=$(echo "$BAM3a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT3e_top=$(echo "$BAM3a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT3e_bottom=$(echo "$BAM3a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT3f_top=$(echo "$BAM3a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT3f_bottom=$(echo "$BAM3a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
SCALE3_OUT=$(echo "$BAM3a" | awk -F. '{print $1"_ForComposite"}')
SCALE3a_OUT=$(echo "$BAM3a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT3_SCALED_COMP_top=$(echo "$BAM3a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
CDT3_SCALED_COMP_bottom=$(echo "$BAM3a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT3_top=$(echo "$BAM3a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
SCALED_OUT3_bottom=$(echo "$BAM3a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
BAM4a=$(echo $BAM4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT4=$(echo "$BAM4a""_""$BEDFILE" | awk -F. '{print $1"_midpoint.out"}')
CDT4=$(echo "$BAM4a""_""$BEDFILE" | awk -F. '{print $1}')
CDT4b=$(echo "$BAM4a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT4c=$(echo "$BAM4a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt"}')
CDT4_SCALED=$(echo "$BAM4a""_""$BEDFILE" | awk -F. '{print $1"_scaled.cdt"}')
SCALE4=$(echo "$BAM4a" | awk -F. '{print $1"_ForCDT"}')
SCALE4a=$(echo "$BAM4a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG4=$(echo "$BAM4a""_""$BEDFILE" | awk -F. '{print $1"_scaled.png"}')
SVG4=$(echo "$BAM4a""_""$BEDFILE" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT4a_top=$(echo "$BAM4a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
OUT4a_bottom=$(echo "$BAM4a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT4d_top=$(echo "$BAM4a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT4d_bottom=$(echo "$BAM4a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT4e_top=$(echo "$BAM4a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT4e_bottom=$(echo "$BAM4a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT4f_top=$(echo "$BAM4a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT4f_bottom=$(echo "$BAM4a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
SCALE4_OUT=$(echo "$BAM4a" | awk -F. '{print $1"_ForComposite"}')
SCALE4a_OUT=$(echo "$BAM4a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT4_SCALED_COMP_top=$(echo "$BAM4a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
CDT4_SCALED_COMP_bottom=$(echo "$BAM4a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT4_top=$(echo "$BAM4a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
SCALED_OUT4_bottom=$(echo "$BAM4a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
BAM5a=$(echo $BAM5 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT5=$(echo "$BAM5a""_""$BEDFILE" | awk -F. '{print $1"_midpoint.out"}')
CDT5=$(echo "$BAM5a""_""$BEDFILE" | awk -F. '{print $1}')
CDT5b=$(echo "$BAM5a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT5c=$(echo "$BAM5a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt"}')
CDT5_SCALED=$(echo "$BAM5a""_""$BEDFILE" | awk -F. '{print $1"_scaled.cdt"}')
SCALE5=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT"}')
SCALE5a=$(echo "$BAM5a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG5=$(echo "$BAM5a""_""$BEDFILE" | awk -F. '{print $1"_scaled.png"}')
SVG5=$(echo "$BAM5a""_""$BEDFILE" | awk -F. '{print $1"_scaled_labeled.svg"}')
OUT5a_top=$(echo "$BAM5a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
OUT5a_bottom=$(echo "$BAM5a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_midpoint.out"}')
CDT5d_top=$(echo "$BAM5a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT5d_bottom=$(echo "$BAM5a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite"}')
CDT5e_top=$(echo "$BAM5a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT5e_bottom=$(echo "$BAM5a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt.gz"}')
CDT5f_top=$(echo "$BAM5a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
CDT5f_bottom=$(echo "$BAM5a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_combined.cdt"}')
SCALE5_OUT=$(echo "$BAM5a" | awk -F. '{print $1"_ForComposite"}')
SCALE5a_OUT=$(echo "$BAM5a" | awk -F. '{print $1"_ForComposite_ScalingFactors.out"}')
CDT5_SCALED_COMP_top=$(echo "$BAM5a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
CDT5_SCALED_COMP_bottom=$(echo "$BAM5a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.cdt"}')
SCALED_OUT5_top=$(echo "$BAM5a""_""$TOP_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')
SCALED_OUT5_bottom=$(echo "$BAM5a""_""$BOTTOM_BED_2200" | awk -F. '{print $1"_ForComposite_scaled.tab"}')

sampleID=Supp_Fig1_v3_240102_ROAR.slurm
rm -f $sampleID
echo "$JOBSTATS" >> $sampleID
echo "#set directory" >> $sampleID
echo "cd $OUTPUT" >> $sampleID
echo "#align bed to reference" >> $sampleID
echo "java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED1 $CpG_BEDFILE $TSS_BEDFILE" >> $sampleID
echo "#sort bedfile by CDT" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -o=K562_CoPRO_TSS_sortByCpGIsland -x=900 2000 $TSS_BEDFILE $ALIGN_BED1" >> $sampleID
echo "#align bed to reference" >> $sampleID
echo "java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED2 $CpG_BEDFILE $BEDFILE" >> $sampleID
echo "#heatmap generation" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$CpG_PNG -c 0000FF -a=1 $ALIGN_BED2" >> $sampleID
echo "#label above heatmaps" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $CpG_PNG --output=$CpG_SVG --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from TSS' --y-label='11,865 CoPRO determined TSSs sorted by CpG island length'" >> $sampleID
echo "#make composite plots for TSS(s) within CpG islands (top10K) and those not in CpG islands (bottom1600)" >> $sampleID
echo "#make bedfiles for top10K sites and bottom 1600 sites" >> $sampleID
echo "cat $BEDFILE | head -2500 > $TOP_BED" >> $sampleID
echo "cat $BEDFILE | tail -1600 > $BOTTOM_BED" >> $sampleID
echo "#expand bedfiles" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2200 $TOP_BED -o=$TOP_BED_2200" >> $sampleID
echo "java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2200 $BOTTOM_BED -o=$BOTTOM_BED_2200" >> $sampleID
echo "#align bed to reference" >> $sampleID
echo "java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED_top $CpG_BEDFILE $TOP_BED_2200" >> $sampleID
echo "java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED_bottom $CpG_BEDFILE $BOTTOM_BED_2200" >> $sampleID
echo "#take average of CDT to make composite plot (now ready for plotter). NOTE: 'sum_Col_CDT.pl' takes the average of each column." >> $sampleID
echo "perl $JOB $ALIGN_BED_top $OUT_top" >> $sampleID
echo "perl $JOB $ALIGN_BED_bottom $OUT_bottom" >> $sampleID
echo "#do analysis for BNase-seq and MNase-seq at TSS bedfile sorted by CpG island length" >> $sampleID
echo "#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $BEDFILE $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2 $BEDFILE $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $BEDFILE $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4 $BEDFILE $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5 $BEDFILE $BAM5" >> $sampleID
echo "#scale output files: options - total tag scaling -t" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE4 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE5 $BAM5" >> $sampleID
echo "#scale BNase-seq, MNase-seq, or DNase-seq data in matrix by scaling factor" >> $sampleID
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
echo "#make heatmaps, merge heatmaps [make benzonase-seq pink and MNase black] **change threshold if necessary" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG1 -c ff00ff -p=0.95 $CDT1_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG2 -p=0.95 $CDT2_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG3 -p=0.95 $CDT3_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG4 -p=0.95 $CDT4_SCALED" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG5 -p=0.95 $CDT5_SCALED" >> $sampleID
echo "#label above heatmaps" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG1 --output=$SVG1 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from TSS (kb)' --y-label='11,865 CoPRO determined TSSs sorted by CpG island length'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG2 --output=$SVG2 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from TSS (kb)' --y-label='11,865 CoPRO determined TSSs sorted by CpG island length'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG3 --output=$SVG3 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from TSS (kb)' --y-label='11,865 CoPRO determined TSSs sorted by CpG island length'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG4 --output=$SVG4 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from TSS (kb)' --y-label='11,865 CoPRO determined TSSs sorted by CpG island length'" >> $sampleID
echo "java -Djava.awt.headless=true -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG5 --output=$SVG5 --width=2 --font-size=18 --left-label='-1' --mid-label=0 --right-label='+1' --x-label='Distance from TSS (kb)' --y-label='11,865 CoPRO determined TSSs sorted by CpG island length'" >> $sampleID
echo "#prep files for plotter" >> $sampleID
echo "#do another tag-pileUp (output is input directory). Settings: midpoint(m), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1d_top -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1a_top $TOP_BED_2200 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1d_bottom -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1a_bottom $BOTTOM_BED_2200 $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2d_top -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2a_top $TOP_BED_2200 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT2d_bottom -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT2a_bottom $BOTTOM_BED_2200 $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3d_top -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3a_top $TOP_BED_2200 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3d_bottom -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3a_bottom $BOTTOM_BED_2200 $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4d_top -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4a_top $TOP_BED_2200 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT4d_bottom -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT4a_bottom $BOTTOM_BED_2200 $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5d_top -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5a_top $TOP_BED_2200 $BAM5" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT5d_bottom -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT5a_bottom $BOTTOM_BED_2200 $BAM5" >> $sampleID
echo "#scale output files: options - total tag scaling -t; *this is same regardless of bedfile" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1_OUT $BAM1" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2_OUT $BAM2" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3_OUT $BAM3" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE4_OUT $BAM4" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE5_OUT $BAM5" >> $sampleID
echo "#scale benzonase-seq and MNase-seq data in matrix by scaling factor" >> $sampleID
echo "gunzip -c $CDT1e_top > $OUTPUT/$CDT1f_top" >> $sampleID
echo "gunzip -c $CDT1e_bottom > $OUTPUT/$CDT1f_bottom" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_COMP_top --scaling-factor=\$(cat $SCALE1a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT1f_top" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED_COMP_bottom --scaling-factor=\$(cat $SCALE1a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT1f_bottom" >> $sampleID
echo "gunzip -c $CDT2e_top > $OUTPUT/$CDT2f_top" >> $sampleID
echo "gunzip -c $CDT2e_bottom > $OUTPUT/$CDT2f_bottom" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_COMP_top --scaling-factor=\$(cat $SCALE2a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT2f_top" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED_COMP_bottom --scaling-factor=\$(cat $SCALE2a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT2f_bottom" >> $sampleID
echo "gunzip -c $CDT3e_top > $OUTPUT/$CDT3f_top" >> $sampleID
echo "gunzip -c $CDT3e_bottom > $OUTPUT/$CDT3f_bottom" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_top --scaling-factor=\$(cat $SCALE3a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT3f_top" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED_COMP_bottom --scaling-factor=\$(cat $SCALE3a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT3f_bottom" >> $sampleID
echo "gunzip -c $CDT4e_top > $OUTPUT/$CDT4f_top" >> $sampleID
echo "gunzip -c $CDT4e_bottom > $OUTPUT/$CDT4f_bottom" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_top --scaling-factor=\$(cat $SCALE4a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT4f_top" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT4_SCALED_COMP_bottom --scaling-factor=\$(cat $SCALE4a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT4f_bottom" >> $sampleID
echo "gunzip -c $CDT5e_top > $OUTPUT/$CDT5f_top" >> $sampleID
echo "gunzip -c $CDT5e_bottom > $OUTPUT/$CDT5f_bottom" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_COMP_top --scaling-factor=\$(cat $SCALE5a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT5f_top" >> $sampleID
echo "java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT5_SCALED_COMP_bottom --scaling-factor=\$(cat $SCALE5a_OUT | cut -f2 | tail -1 | awk '{print \$1}') $CDT5f_bottom" >> $sampleID
echo "#make scale OUT file" >> $sampleID
echo "perl $JOB $CDT1_SCALED_COMP_top $SCALED_OUT1_top" >> $sampleID
echo "perl $JOB $CDT1_SCALED_COMP_bottom $SCALED_OUT1_bottom" >> $sampleID
echo "perl $JOB $CDT2_SCALED_COMP_top $SCALED_OUT2_top" >> $sampleID
echo "perl $JOB $CDT2_SCALED_COMP_bottom $SCALED_OUT2_bottom" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_top $SCALED_OUT3_top" >> $sampleID
echo "perl $JOB $CDT3_SCALED_COMP_bottom $SCALED_OUT3_bottom" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_top $SCALED_OUT4_top" >> $sampleID
echo "perl $JOB $CDT4_SCALED_COMP_bottom $SCALED_OUT4_bottom" >> $sampleID
echo "perl $JOB $CDT5_SCALED_COMP_top $SCALED_OUT5_top" >> $sampleID
echo "perl $JOB $CDT5_SCALED_COMP_bottom $SCALED_OUT5_bottom" >> $sampleID
echo "#remove intermediate files" >> $sampleID
echo "rm $TSS_BEDFILEa" >> $sampleID
echo "rm $CpG_BEDFILEa" >> $sampleID
echo "rm $BEDFILE" >> $sampleID
echo "rm $BEDFILEb" >> $sampleID
echo "rm $ALIGN_BED1" >> $sampleID
echo "rm $ALIGN_BED2" >> $sampleID
echo "rm $CpG_PNG" >> $sampleID
echo "rm $TOP_BED" >> $sampleID
echo "rm $BOTTOM_BED" >> $sampleID
echo "rm $TOP_BED_2200" >> $sampleID
echo "rm $BOTTOM_BED_2200" >> $sampleID
echo "rm $ALIGN_BED_top" >> $sampleID
echo "rm $ALIGN_BED_bottom" >> $sampleID
echo "rm $OUT1" >> $sampleID
echo "rm $CDT1b" >> $sampleID
echo "rm $CDT1c" >> $sampleID
echo "rm $CDT1_SCALED" >> $sampleID
echo "rm $SCALE1a" >> $sampleID
echo "rm $PNG1" >> $sampleID
echo "rm $OUT1a_top" >> $sampleID
echo "rm $OUT1a_bottom" >> $sampleID
echo "rm $CDT1e_top" >> $sampleID
echo "rm $CDT1e_bottom" >> $sampleID
echo "rm $CDT1f_top" >> $sampleID
echo "rm $CDT1f_bottom" >> $sampleID
echo "rm $SCALE1a_OUT" >> $sampleID
echo "rm $CDT1_SCALED_COMP_top" >> $sampleID
echo "rm $CDT1_SCALED_COMP_bottom" >> $sampleID
echo "rm $OUT2" >> $sampleID
echo "rm $CDT2b" >> $sampleID
echo "rm $CDT2c" >> $sampleID
echo "rm $CDT2_SCALED" >> $sampleID
echo "rm $SCALE2a" >> $sampleID
echo "rm $PNG2" >> $sampleID
echo "rm $OUT2a_top" >> $sampleID
echo "rm $OUT2a_bottom" >> $sampleID
echo "rm $CDT2e_top" >> $sampleID
echo "rm $CDT2e_bottom" >> $sampleID
echo "rm $CDT2f_top" >> $sampleID
echo "rm $CDT2f_bottom" >> $sampleID
echo "rm $SCALE2a_OUT" >> $sampleID
echo "rm $CDT2_SCALED_COMP_top" >> $sampleID
echo "rm $CDT2_SCALED_COMP_bottom" >> $sampleID
echo "rm $OUT3" >> $sampleID
echo "rm $CDT3b" >> $sampleID
echo "rm $CDT3c" >> $sampleID
echo "rm $CDT3_SCALED" >> $sampleID
echo "rm $SCALE3a" >> $sampleID
echo "rm $PNG3" >> $sampleID
echo "rm $OUT3a_top" >> $sampleID
echo "rm $OUT3a_bottom" >> $sampleID
echo "rm $CDT3e_top" >> $sampleID
echo "rm $CDT3e_bottom" >> $sampleID
echo "rm $CDT3f_top" >> $sampleID
echo "rm $CDT3f_bottom" >> $sampleID
echo "rm $SCALE3a_OUT" >> $sampleID
echo "rm $CDT3_SCALED_COMP_top" >> $sampleID
echo "rm $CDT3_SCALED_COMP_bottom" >> $sampleID
echo "rm $OUT4" >> $sampleID
echo "rm $CDT4b" >> $sampleID
echo "rm $CDT4c" >> $sampleID
echo "rm $CDT4_SCALED" >> $sampleID
echo "rm $SCALE4a" >> $sampleID
echo "rm $PNG4" >> $sampleID
echo "rm $OUT4a_top" >> $sampleID
echo "rm $OUT4a_bottom" >> $sampleID
echo "rm $CDT4e_top" >> $sampleID
echo "rm $CDT4e_bottom" >> $sampleID
echo "rm $CDT4f_top" >> $sampleID
echo "rm $CDT4f_bottom" >> $sampleID
echo "rm $SCALE4a_OUT" >> $sampleID
echo "rm $CDT4_SCALED_COMP_top" >> $sampleID
echo "rm $CDT4_SCALED_COMP_bottom" >> $sampleID
echo "rm $OUT5" >> $sampleID
echo "rm $CDT5b" >> $sampleID
echo "rm $CDT5c" >> $sampleID
echo "rm $CDT5_SCALED" >> $sampleID
echo "rm $SCALE5a" >> $sampleID
echo "rm $PNG5" >> $sampleID
echo "rm $OUT5a_top" >> $sampleID
echo "rm $OUT5a_bottom" >> $sampleID
echo "rm $CDT5e_top" >> $sampleID
echo "rm $CDT5e_bottom" >> $sampleID
echo "rm $CDT5f_top" >> $sampleID
echo "rm $CDT5f_bottom" >> $sampleID
echo "rm $SCALE5a_OUT" >> $sampleID
echo "rm $CDT5_SCALED_COMP_top" >> $sampleID
echo "rm $CDT5_SCALED_COMP_bottom" >> $sampleID
echo "## finish script" >> $sampleID
