# purpose - make Fig. 1 showing benzonase-seq and MNase-seq at ref-seq TSS(s) sorted by K562 gene expression.

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
TSS_BEDFILE=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig1_atTSS_CpGsort/bedfiles/K562_CoPRO-expressed_Gene-refSeqTSS_2000bp.bed
CpG_BEDFILE=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig1_atTSS_CpGsort/bedfiles/UCSCgb_hg19_CpGislands_230426.bed

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig1_atTSS_CpGsort/Fig1_230926_output

#set bam library file to BI_rep1
BAM1=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/K562_benzonase-seq_master.bam
BAM2=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/00_REFERENCE_DATA/K562_Datasets/MNase/K562_MNase.bam
BAM3=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files/SRR16815400_master.bam

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
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
OUT2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_read1.out"}')
CDT2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1}')
CDT2b=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt.gz"}')
CDT2c=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_combined.cdt"}')
CDT2_SCALED=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_scaled.cdt"}')
SCALE2=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT"}')
SCALE2a=$(echo "$BAM2a" | awk -F. '{print $1"_ForCDT_ScalingFactors.out"}')
PNG2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_scaled.png"}')
SVG2=$(echo "$BAM2a""_""$BEDFILE" | awk -F. '{print $1"_scaled_labeled.svg"}')
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

#set directory
cd $OUTPUT

#align bed to reference
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED1 $CpG_BEDFILE $TSS_BEDFILE

#sort bedfile by CDT
java -jar $SCRIPTMANAGER coordinate-manipulation sort-bed -o=K562_CoPRO_TSS_sortByCpGIsland -x=900 2000 $TSS_BEDFILE $ALIGN_BED1

#align bed to reference
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED2 $CpG_BEDFILE $BEDFILE

#heatmap generation
java -jar $SCRIPTMANAGER figure-generation heatmap -o=$CpG_PNG -c 0000FF -a=1 $ALIGN_BED2

#label above heatmaps
java -jar $SCRIPTMANAGER figure-generation label-heatmap $CpG_PNG --output=$CpG_SVG --width=2 --font-size=18 --left-label="-1" --mid-label=0 --right-label="+1" --x-label="Distance from TSS (kb)" --y-label="11,865 CoPRO determined TSSs sorted by CpG island length"

#make composite plots for TSS(s) within CpG islands (top10K) and those not in CpG islands (bottom1600)
#make bedfiles for top10K sites and bottom 1600 sites
cat $BEDFILE | head -2500 > $TOP_BED
cat $BEDFILE | tail -1600 > $BOTTOM_BED

#expand bedfiles
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2200 $TOP_BED -o=$TOP_BED_2200
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c=2200 $BOTTOM_BED -o=$BOTTOM_BED_2200

#align bed to reference
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED_top $CpG_BEDFILE $TOP_BED_2200
java -jar $SCRIPTMANAGER peak-analysis peak-align-ref -o=$ALIGN_BED_bottom $CpG_BEDFILE $BOTTOM_BED_2200

#take average of CDT to make composite plot (now ready for plotter). NOTE: 'sum_Col_CDT.pl' takes the average of each column.
perl $JOB $ALIGN_BED_top $OUT_top
perl $JOB $ALIGN_BED_bottom $OUT_bottom

#do analysis for BNase-seq and MNase-seq at TSS bedfile sorted by CpG island length
#do initial tag-pileUp (output is input directory). Settings: midpoint(m) OR 5 prime end (-5) with read 1 (-1), Gizp output cdt (z), No smoothing (N), required proper PEs (p), load blacklist **total tag option (-t) removed**
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT1 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT1 $BEDFILE $BAM1
java -jar $SCRIPTMANAGER read-analysis tag-pileup -5 -1 -z --combined --output-matrix=$CDT2 -N --cpu=4 --shift=80 --blacklist-filter=$BLACKLIST -o=$OUT2 $BEDFILE $BAM2
java -jar $SCRIPTMANAGER read-analysis tag-pileup -m -z --combined --output-matrix=$CDT3 -N -p --cpu=4 --blacklist-filter=$BLACKLIST -o=$OUT3 $BEDFILE $BAM3

#scale output files: options - total tag scaling -t
java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE1 $BAM1
java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE2 $BAM2
java -jar $SCRIPTMANAGER read-analysis scaling-factor -t --blacklist=$BLACKLIST -o=$SCALE3 $BAM3

#scale benzonase-seq and MNase-seq data in matrix by scaling factor
gunzip -c $CDT1b > $OUTPUT/$CDT1c
SCALE1b=$(cat $SCALE1a | cut -f 2 | tail -1 | awk '{print $1}')
java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT1_SCALED --scaling-factor=$SCALE1b $CDT1c
gunzip -c $CDT2b > $OUTPUT/$CDT2c
SCALE2b=$(cat $SCALE2a | cut -f 2 | tail -1 | awk '{print $1}')
java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT2_SCALED --scaling-factor=$SCALE2b $CDT2c
gunzip -c $CDT3b > $OUTPUT/$CDT3c
SCALE3b=$(cat $SCALE3a | cut -f 2 | tail -1 | awk '{print $1}')
java -jar $SCRIPTMANAGER read-analysis scale-matrix -o=$CDT3_SCALED --scaling-factor=$SCALE3b $CDT3c

#make heatmaps, merge heatmaps [make benzonase-seq pink and MNase black] **change threshold if necessary
java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG1 -c ff00ff -p=0.95 $CDT1_SCALED
java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG2 -p=0.95 $CDT2_SCALED
java -jar $SCRIPTMANAGER figure-generation heatmap -o=$PNG3 -c 00ff1e -p=0.95 $CDT3_SCALED

#label above heatmaps
java -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG1 --output=$SVG1 --width=2 --font-size=18 --left-label="-1" --mid-label=0 --right-label="+1" --x-label="Distance from TSS (kb)" --y-label="11,865 CoPRO determined TSSs sorted by CpG island length"
java -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG2 --output=$SVG2 --width=2 --font-size=18 --left-label="-1" --mid-label=0 --right-label="+1" --x-label="Distance from TSS (kb)" --y-label="11,865 CoPRO determined TSSs sorted by CpG island length"
java -jar $SCRIPTMANAGER figure-generation label-heatmap $PNG3 --output=$SVG3 --width=2 --font-size=18 --left-label="-1" --mid-label=0 --right-label="+1" --x-label="Distance from TSS (kb)" --y-label="11,865 CoPRO determined TSSs sorted by CpG island length"

# finish script
echo "DONE"
