# purpose - make Fig. 4c showing PE insert size histogram of BNase-seq.

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

#output directory
OUTPUT=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig4_InsertSizeHistogram/histogram

#set bam library file to BI_rep1
BAM1=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230718_MERGE/K562_benzonase-seq_master.bam

#set scriptmanager and job
SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig4_InsertSizeHistogram/jobs/make_fragment_histograms-MOD.py
JOB2=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig4_InsertSizeHistogram/jobs/make_fragment_histograms-MOD2.py

#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

#set output file names
BAM1a=$(echo $BAM1 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HISTOGRAM1=$(echo "$BAM1a" | awk -F. '{print $1}')
HISTOGRAM1b=$(echo "$BAM1a" | awk -F. '{print $1"_InsertHistogram.out"}')
COMP1=$(echo "$BAM1a" | awk -F. '{print $1"_compositePlot.svg"}')
COMP2=$(echo "$BAM1a" | awk -F. '{print $1"_compositePlot_20bp.svg"}')

#set directory
cd $OUTPUT

#expand bedfile in scriptmanager (output is input directory) **add variable
java -jar $SCRIPTMANAGER bam-statistics pe-stat $BAM1 --min=0 -o=$HISTOGRAM1 --max=340

#make frequency histograms
python $JOB -i $HISTOGRAM1b -o $COMP1
python $JOB2 -i $HISTOGRAM1b -o $COMP2

# finish script
echo "DONE"
