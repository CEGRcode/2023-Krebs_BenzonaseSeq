# purpose - make Fig. 4a and Fig. 4b showing PE insert size histograms of DNase-seq and MNase-seq libraries.

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
BAM2=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files/SRR16815400_master.bam
BAM3=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files/SRR3211679_master.bam
BAM4=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files/SRR3211680_master.bam
BAM5=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files/SRR3211681_master.bam
BAM6=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/230810_MNase_DNase/final_files/SRR3211682_master.bam

#set scriptmanager and job
SCRIPTMANAGER=/gpfs/group/bfp2/default/pughlab-members/juk398-JordanKrebs/scriptmanager/build/libs/ScriptManager-v0.14.jar
JOB=/gpfs/group/bfp2/default/pughlab-members/wkl2-WillLai/NucleosomeAtlas_Project/figures/fig4_InsertSizeHistogram/jobs/make_fragment_histograms-MOD.py


#------ CODE ------

# stop on errors & undefined variables, print commands
# defense against the dark arts
set -eux
echo "defense against the dark arts activated"

#set output file names
BAM2a=$(echo $BAM2 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HISTOGRAM2=$(echo "$BAM2a" | awk -F. '{print $1}')
HISTOGRAM2b=$(echo "$BAM2a" | awk -F. '{print $1"_InsertHistogram.out"}')
COMP2=$(echo "$BAM2a" | awk -F. '{print $1"_compositePlot.svg"}')
BAM3a=$(echo $BAM3 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HISTOGRAM3=$(echo "$BAM3a" | awk -F. '{print $1}')
HISTOGRAM3b=$(echo "$BAM3a" | awk -F. '{print $1"_InsertHistogram.out"}')
COMP3=$(echo "$BAM3a" | awk -F. '{print $1"_compositePlot.svg"}')
BAM4a=$(echo $BAM4 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HISTOGRAM4=$(echo "$BAM4a" | awk -F. '{print $1}')
HISTOGRAM4b=$(echo "$BAM4a" | awk -F. '{print $1"_InsertHistogram.out"}')
COMP4=$(echo "$BAM4a" | awk -F. '{print $1"_compositePlot.svg"}')
BAM5a=$(echo $BAM5 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HISTOGRAM5=$(echo "$BAM5a" | awk -F. '{print $1}')
HISTOGRAM5b=$(echo "$BAM5a" | awk -F. '{print $1"_InsertHistogram.out"}')
COMP5=$(echo "$BAM5a" | awk -F. '{print $1"_compositePlot.svg"}')
BAM6a=$(echo $BAM6 | rev | cut -d"/" -f1 | rev | awk -F. '{print $1}')
HISTOGRAM6=$(echo "$BAM6a" | awk -F. '{print $1}')
HISTOGRAM6b=$(echo "$BAM6a" | awk -F. '{print $1"_InsertHistogram.out"}')
COMP6=$(echo "$BAM6a" | awk -F. '{print $1"_compositePlot.svg"}')

#set directory
cd $OUTPUT

#expand bedfile in scriptmanager (output is input directory) **add variable
java -jar $SCRIPTMANAGER bam-statistics pe-stat $BAM2 --min=0 -o=$HISTOGRAM2 --max=340
java -jar $SCRIPTMANAGER bam-statistics pe-stat $BAM3 --min=0 -o=$HISTOGRAM3 --max=340
java -jar $SCRIPTMANAGER bam-statistics pe-stat $BAM4 --min=0 -o=$HISTOGRAM4 --max=340
java -jar $SCRIPTMANAGER bam-statistics pe-stat $BAM5 --min=0 -o=$HISTOGRAM5 --max=340
java -jar $SCRIPTMANAGER bam-statistics pe-stat $BAM6 --min=0 -o=$HISTOGRAM6 --max=340

#make frequency histograms
python $JOB -i $HISTOGRAM2b -o $COMP2
python $JOB -i $HISTOGRAM3b -o $COMP3
python $JOB -i $HISTOGRAM4b -o $COMP4
python $JOB -i $HISTOGRAM5b -o $COMP5
python $JOB -i $HISTOGRAM6b -o $COMP6

# finish script
echo "DONE"
